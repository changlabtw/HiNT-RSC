"""Command-line interface: hint-rsc {run, sweep, coverage, covariates, evaluate, check}."""
import argparse
import shutil
import subprocess
import sys

from . import __version__
from .smooth import METHODS

EPILOG = '''examples:
  hint-rsc run --hic sample.hic --genome hg38 --enzyme HindIII -o out/            # HiNT-CNV baseline, lambda 3
  hint-rsc sweep --coverage out/coverage.tsv --genome hg38 --enzyme HindIII -o out/sweep   # tune lambda (and w)
  hint-rsc run --coverage out/coverage.tsv --genome hg38 --enzyme HindIII --lambda 1 -o out/l1
  hint-rsc covariates --fasta hg38.fa.gz --genome hg38 --enzyme DpnII --mappability k50.Umap.bedGraph.gz -o hg38_DpnII.tsv

Recommendations from the paper: tune lambda per dataset (sweep), exclude hg38 problem regions (default on hg38), and
smooth when a bin-level copy-number profile is the goal (profile.tsv is smoothed with --profile-window, default 9).
'''


def _pipeline_args(p, sweep=False):
    g = p.add_argument_group('input')
    g.add_argument('--hic', help='.hic file (coverage is computed and cached as <outdir>/coverage.tsv)')
    g.add_argument('--coverage', help='coverage table from "hint-rsc coverage" (instead of --hic)')
    g.add_argument('--genome', required=True, help='hg19, hg38, or a chromosome-sizes file (chrom<TAB>length)')
    g.add_argument('--enzyme', help='restriction enzyme, selects the bundled covariates (hg19 MboI, hg38 HindIII)')
    g.add_argument('--covariates', help='covariate table from "hint-rsc covariates" (overrides --enzyme)')
    g.add_argument('--contacts', choices=('complete', 'intra'), default='complete',
                   help='complete (intra + inter) or intra-chromosomal contacts only (faster; default: %(default)s)')
    g.add_argument('--resolution', type=int, default=50000, help='bin size in bp (default: %(default)s; the paper used 50 kb)')
    g.add_argument('-o', '--outdir', required=True)
    g = p.add_argument_group('normalization')
    g.add_argument('--problem-regions', nargs='+', default=['auto'], metavar='FILE',
                   help="'auto' (bundled hg38 centromeres, gaps and ENCODE blacklist v2; nothing for other genomes), 'none', "
                        'or BED / UCSC cytoBand / UCSC gap files (default: auto)')
    g.add_argument('--min-frac', type=float, default=0.25, help='exclude a bin when this fraction lies in problem regions (default: %(default)s)')
    g.add_argument('--control-hic', help='control .hic file for control correction (residual - f x control residual)')
    g.add_argument('--control-coverage', help='control coverage table (instead of --control-hic)')
    g.add_argument('--control-enzyme', help='control enzyme if it differs from --enzyme')
    g.add_argument('--control-covariates', help='control covariate table if it differs from --covariates')
    g.add_argument('--control-fraction', type=float, default=0.3, help='subtraction fraction f (default: %(default)s)')
    g = p.add_argument_group('smoothing and segmentation')
    g.add_argument('--smoother', choices=METHODS, default='mean', help='smoothing filter (default: %(default)s)')
    if sweep:
        g.add_argument('--windows', default='0,3,5,9,21', help='comma-separated windows; 0 = no smoothing (default: %(default)s)')
        g.add_argument('--lambdas', default='0.5,1,1.5,2,3,4,5', help='comma-separated BIC-seq2 penalties (default: %(default)s)')
        g.add_argument('--workers', type=int, default=2, help='parallel BIC-seq2 runs (default: %(default)s)')
    else:
        g.add_argument('-w', '--window', type=int, default=0,
                       help='smoothing window in bins before segmentation; 0 = none, as in HiNT-CNV (default: %(default)s)')
        g.add_argument('--lambda', dest='lam', type=float, default=3,
                       help='BIC-seq2 penalty (default: %(default)s, the HiNT-CNV default; tune it with "hint-rsc sweep")')
        g.add_argument('--profile-window', type=int, default=9,
                       help='mean-filter window for the bin-level profile.tsv; 0 = none (default: %(default)s)')
    g.add_argument('--threshold', type=float, default=0.3, help='CNV call threshold on |log2 ratio| (default: %(default)s)')
    g = p.add_argument_group('reference (optional)')
    g.add_argument('--reference', help='reference copy-number segments to score against (e.g. from WGS)')
    g.add_argument('--reference-format', choices=('auto', 'seg', 'bed'), default='auto',
                   help='seg: BIC-seq2-style with header, 1-based; bed: chrom start end log2 (default: auto by extension)')
    g = p.add_argument_group('tools')
    g.add_argument('--bicseq', help='NBICseq-seg.pl (default: $BICSEQ, then PATH)')
    g.add_argument('--rscript', default='Rscript', help='Rscript with mgcv and optparse (default: %(default)s)')


def cmd_check(args):
    ok = True
    for mod in ('numpy', 'scipy', 'pandas', 'pywt', 'hicstraw'):
        try:
            __import__(mod)
            print(f'ok      python module {mod}')
        except ImportError:
            ok = False
            print(f'MISSING python module {mod}' + (' (needed only for .hic input)' if mod == 'hicstraw' else ''))
    r = shutil.which(args.rscript)
    if r and subprocess.run([r, '-e', 'library(mgcv); library(optparse)'], capture_output=True).returncode == 0:
        print(f'ok      {r} with mgcv and optparse')
    else:
        ok = False
        print(f'MISSING Rscript with mgcv and optparse ({args.rscript})')
    from .segment import find_bicseq
    try:
        print(f'ok      BIC-seq2 {find_bicseq(args.bicseq)}')
    except SystemExit as e:
        ok = False
        print(f'MISSING {e}')
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(prog='hint-rsc', description='HiNT-RSC: residual smoothing for Hi-C-based CNV detection '
                                 '(HiNT-CNV GAM normalization + BIC-seq2 segmentation).', epilog=EPILOG,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version', action='version', version=f'hint-rsc {__version__}')
    sub = ap.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('run', help='one setting: CNV segments, calls and a bin-level profile',
                       description='GAM residuals -> optional control correction and smoothing -> BIC-seq2. Writes segments.seg, '
                                   'cnv_calls.bed, profile.tsv, run.json (and metrics.tsv with --reference) to --outdir.')
    _pipeline_args(p)
    p = sub.add_parser('sweep', help='grid of smoothing windows x lambdas, to tune the segmentation granularity',
                       description='Segments every window x lambda combination and writes sweep.tsv (calls, called Mb and, '
                                   'with --reference, the agreement metrics) to --outdir.')
    _pipeline_args(p, sweep=True)

    p = sub.add_parser('coverage', help='1D coverage (row sums) from a .hic file')
    p.add_argument('--hic', required=True)
    p.add_argument('--genome', required=True, help='hg19, hg38, or a chromosome-sizes file')
    p.add_argument('--resolution', type=int, default=50000)
    p.add_argument('--intra-only', action='store_true', help='intra-chromosomal contacts only (faster)')
    p.add_argument('--chroms', help='comma-separated subset of chromosomes')
    p.add_argument('-o', '--output', required=True)

    p = sub.add_parser('covariates', help='GAM covariates (restriction sites, GC, mappability) from a genome FASTA')
    p.add_argument('--fasta', required=True, help='genome FASTA (.fa or .fa.gz) matching --genome')
    p.add_argument('--genome', required=True, help='hg19, hg38, or a chromosome-sizes file')
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--enzyme', help='HindIII, MboI, DpnII, NcoI or Arima (GATC + GANTC)')
    g.add_argument('--motif', help='restriction motif(s), comma-separated, N = any base (e.g. GATC,GANTC)')
    p.add_argument('--mappability', help='mappability bedGraph (.gz), e.g. Umap k50 converted with bigWigToBedGraph; omit for 0')
    p.add_argument('--resolution', type=int, default=50000)
    p.add_argument('-o', '--output', required=True)

    p = sub.add_parser('evaluate', help='score a segmentation against a reference copy-number profile')
    p.add_argument('--seg', required=True, help='Hi-C segmentation (segments.seg or a sweep segmentation)')
    p.add_argument('--bins', help='its BIC-seq2 bin directory, for the bin-level correlations')
    p.add_argument('--reference', required=True)
    p.add_argument('--reference-format', choices=('auto', 'seg', 'bed'), default='auto')
    p.add_argument('--genome', required=True)
    p.add_argument('--resolution', type=int, default=50000)
    p.add_argument('--threshold', type=float, default=0.3)
    p.add_argument('--iou', type=float, default=0.5, help='matching threshold on intersection over union (default: %(default)s)')

    p = sub.add_parser('check', help='check that R (mgcv), BIC-seq2 and the Python modules are available')
    p.add_argument('--bicseq')
    p.add_argument('--rscript', default='Rscript')

    args = ap.parse_args(argv)
    if args.cmd == 'run':
        from .pipeline import run
        return run(args)
    if args.cmd == 'sweep':
        from .pipeline import sweep
        return sweep(args)
    if args.cmd == 'check':
        return cmd_check(args)
    from .genome import Genome
    genome = Genome(args.genome, args.resolution)
    if args.cmd == 'coverage':
        from .coverage import coverage_profiles, write_coverage
        chroms = args.chroms.split(',') if args.chroms else None
        log = lambda m: print(m, file=sys.stderr, flush=True)  # noqa: E731
        write_coverage(coverage_profiles(args.hic, genome, chroms, args.intra_only, log), args.output)
    elif args.cmd == 'covariates':
        from .covariates import ENZYME_MOTIFS, build_covariates
        if args.enzyme and args.enzyme not in ENZYME_MOTIFS:
            raise SystemExit(f'unknown enzyme {args.enzyme}; use --motif (known: {", ".join(ENZYME_MOTIFS)})')
        build_covariates(args.fasta, genome, args.motif or ENZYME_MOTIFS[args.enzyme], args.mappability, args.output,
                         log=lambda m: print(m, file=sys.stderr, flush=True))
    elif args.cmd == 'evaluate':
        from .evaluate import evaluate, read_reference
        row = evaluate(args.seg, args.bins, read_reference(args.reference, genome, args.reference_format), genome,
                       args.threshold, args.iou)
        print('\t'.join(row))
        print('\t'.join(str(v) for v in row.values()))
    return 0


if __name__ == '__main__':
    sys.exit(main())
