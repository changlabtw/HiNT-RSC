"""HiNT-RSC pipeline: coverage -> GAM residuals -> (control correction) -> smoothing -> min-shift bins -> BIC-seq2."""
import json
import os
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from . import __version__
from .coverage import coverage_profiles, write_coverage
from .evaluate import evaluate, read_hic_seg, read_reference
from .genome import Genome
from .normalize import exclusion_mask, fit_gam, load_table, read_problem_regions
from .segment import bicseq, find_bicseq, fmt, write_bins, write_calls, write_profile
from .smooth import smooth_dict


def log(msg):
    print(f'[hint-rsc] {msg}', file=sys.stderr, flush=True)


def _coverage_file(args, genome, hic, out_name, contacts):
    """Coverage table for args (computed from a .hic file into the output folder when needed)."""
    need = 'rowsum_complete' if contacts == 'complete' else 'rowsum_intra'
    out = os.path.join(args.outdir, out_name)
    if os.path.exists(out) and need in open(out).readline().split():
        log(f'using cached coverage {out}')
        return out
    log(f'computing {"intra-chromosomal" if contacts == "intra" else "complete"} coverage from {hic}')
    return write_coverage(coverage_profiles(hic, genome, intra_only=contacts == 'intra', log=log), out)


def _covariates(spec, genome, enzyme):
    if spec:
        return spec
    path = genome.builtin_covariates(enzyme) if enzyme else None
    if path and genome.resolution != 50000:
        raise SystemExit('the bundled covariates are for 50-kb bins; build a table for this resolution with "hint-rsc covariates"')
    if path is None:
        raise SystemExit(f'no bundled covariates for {genome.name} {enzyme or "(no --enzyme given)"}: bundled are hg19 MboI '
                         'and hg38 HindIII. Build a table with "hint-rsc covariates" and pass it with --covariates.')
    return path


def _exclusion(args, genome):
    if args.problem_regions == ['none']:
        return None, 'none'
    if args.problem_regions == ['auto']:
        files = genome.builtin_problem_regions()
        if files is None:
            log(f'no bundled problem regions for {genome.name}: no bins are excluded (pass --problem-regions BED to exclude)')
            return None, 'none'
        source = f'bundled {genome.name} (cytoBand acen, gap, ENCODE blacklist v2)'
    else:
        files, source = args.problem_regions, ', '.join(args.problem_regions)
    mask = exclusion_mask(read_problem_regions(files), genome, args.min_frac)
    log(f'excluding {sum(int(m.sum()) for m in mask.values())} bins (>= {args.min_frac:g} in problem regions: {source})')
    return mask, source


def prepare(args):
    """Fit the GAM (and the control GAM) and return the genome, the residuals, the fitted mask and run metadata."""
    genome = Genome(args.genome, args.resolution)
    os.makedirs(args.outdir, exist_ok=True)
    if bool(args.hic) == bool(args.coverage):
        raise SystemExit('give one of --hic or --coverage')
    coverage = args.coverage or _coverage_file(args, genome, args.hic, 'coverage.tsv', args.contacts)
    covariates = _covariates(args.covariates, genome, args.enzyme)
    excluded, excl_source = _exclusion(args, genome)
    column = 'rowsum_complete' if args.contacts == 'complete' else 'rowsum_intra'
    resid, fitted = fit_gam(load_table(coverage, covariates), column, genome, os.path.join(args.outdir, 'gam', args.contacts),
                            args.rscript, excluded)
    log(f'GAM fitted on {sum(int(m.sum()) for m in fitted.values())} bins ({args.contacts} contacts)')
    meta = {'version': __version__, 'genome': genome.name, 'resolution': genome.resolution, 'contacts': args.contacts,
            'coverage': coverage, 'covariates': covariates, 'problem_regions': excl_source, 'min_frac': args.min_frac}
    if args.control_hic or args.control_coverage:
        ccov = args.control_coverage or _coverage_file(args, genome, args.control_hic, 'control_coverage.tsv', args.contacts)
        ccovariates = _covariates(args.control_covariates, genome, args.control_enzyme or args.enzyme) \
            if (args.control_covariates or args.control_enzyme) else covariates
        cres, _ = fit_gam(load_table(ccov, ccovariates), column, genome, os.path.join(args.outdir, 'gam', 'control_' + args.contacts),
                          args.rscript, excluded)
        f = args.control_fraction
        resid = {c: resid[c] - f * cres[c] for c in genome.nbins}
        log(f'control correction: residual - {f:g} x control residual')
        meta.update({'control_coverage': ccov, 'control_covariates': ccovariates, 'control_fraction': f})
    return genome, resid, fitted, meta


def _tag(args, window):
    tag = 'none_0' if not window else f'{args.smoother}_{window}'
    return (f'ctrl{fmt(args.control_fraction)}_' + tag) if 'control_fraction' in args.meta else tag


def segment_one(args, genome, resid, window, lam, bicseq_pl):
    bin_dir = os.path.join(args.outdir, 'bins', _tag(args, window))
    config = write_bins(smooth_dict(resid, args.smoother, window), bin_dir, genome.resolution)
    seg = bicseq(config, lam, bicseq_pl, os.path.join(args.outdir, 'tmp', f'{os.path.basename(bin_dir)}_l{fmt(lam)}'))
    return bin_dir, seg


def run(args):
    bicseq_pl = find_bicseq(args.bicseq)
    genome, resid, fitted, args.meta = prepare(args)
    bin_dir, seg = segment_one(args, genome, resid, args.window, args.lam, bicseq_pl)
    shutil.copyfile(seg, os.path.join(args.outdir, 'segments.seg'))
    n_calls = write_calls(seg, genome, os.path.join(args.outdir, 'cnv_calls.bed'), args.threshold)
    write_profile(smooth_dict(resid, 'mean', args.profile_window), fitted, genome, os.path.join(args.outdir, 'profile.tsv'))
    meta = dict(args.meta, smoother=args.smoother, window=args.window, profile_mean_window=args.profile_window, lam=args.lam,
                threshold=args.threshold, bins=bin_dir, segments=len(read_hic_seg(seg, genome)), cnv_calls=n_calls)
    if args.reference:
        row = evaluate(seg, bin_dir, read_reference(args.reference, genome, args.reference_format), genome, args.threshold)
        pd.DataFrame([row]).to_csv(os.path.join(args.outdir, 'metrics.tsv'), sep='\t', index=False)
        meta['metrics'] = row
        log('metrics: ' + ', '.join(f'{k} {v}' for k, v in row.items()))
    with open(os.path.join(args.outdir, 'run.json'), 'w') as f:
        json.dump(meta, f, indent=2)
    log(f'{meta["segments"]} segments, {n_calls} CNV calls (|log2| >= {args.threshold:g}, lambda {fmt(args.lam)}, '
        f'{"no smoothing" if not args.window else f"{args.smoother} w = {args.window}"}) -> {args.outdir}')


def sweep(args):
    bicseq_pl = find_bicseq(args.bicseq)
    genome, resid, _, args.meta = prepare(args)
    windows = [int(w) for w in args.windows.split(',')]
    lambdas = [float(l) for l in args.lambdas.split(',')]
    reference = read_reference(args.reference, genome, args.reference_format) if args.reference else None
    for w in windows:                                     # write each bin set once, before the parallel segmentation
        write_bins(smooth_dict(resid, args.smoother, w), os.path.join(args.outdir, 'bins', _tag(args, w)), genome.resolution)
    jobs = [(w, l) for w in windows for l in lambdas]
    log(f'{len(jobs)} segmentations ({len(windows)} windows x {len(lambdas)} lambdas, {args.workers} workers)')

    def one(job):
        w, lam = job
        bin_dir = os.path.join(args.outdir, 'bins', _tag(args, w))
        seg = bicseq(os.path.join(bin_dir, '_config.txt'), lam, bicseq_pl, os.path.join(args.outdir, 'tmp', f'{_tag(args, w)}_l{fmt(lam)}'))
        segs = read_hic_seg(seg, genome)
        called = [s for s in segs if abs(s[3]) >= args.threshold]
        row = {'smoother': args.smoother if w else 'none', 'window': w, 'lambda': fmt(lam), 'segments': len(segs), 'calls': len(called),
               'gains': sum(s[3] > 0 for s in called), 'losses': sum(s[3] < 0 for s in called),
               'called_Mb': round(sum(s[2] - s[1] + 1 for s in called) / 1e6, 1)}
        if reference is not None:
            row.update(evaluate(seg, bin_dir, reference, genome, args.threshold))
        row['segmentation'] = os.path.relpath(seg, args.outdir)
        return row
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        rows = list(ex.map(one, jobs))
    df = pd.DataFrame(rows)
    out = os.path.join(args.outdir, 'sweep.tsv')
    df.to_csv(out, sep='\t', index=False)
    with open(os.path.join(args.outdir, 'run.json'), 'w') as f:
        json.dump(dict(args.meta, smoother=args.smoother, windows=windows, lambdas=lambdas, threshold=args.threshold,
                       reference=args.reference), f, indent=2)
    print(df.drop(columns=['segmentation', 'matched_by_ref_size'], errors='ignore').to_string(index=False))
    log(f'wrote {out}')
