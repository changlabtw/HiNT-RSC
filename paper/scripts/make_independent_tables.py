"""Tables for the independent cell lines (PANC-1, Caki2): main-text Table 7 and Supplementary Tables S16, S17 and S18.

Inputs (run from paper/; defaults are the outputs of workflow/20_run_independent.sh, or expected/ with --expected):
  PANC-1 / Caki2 with problem regions excluded   panc1_rerun.py --line panc1|caki2       -> output/<line>/excl/results.tsv
  PANC-1 without exclusion (pre-registered)      run_validation.py + score_reference.py  -> output/panc1/noexcl/
Writes <outdir>/Table7_independent_lines.tsv, TableS16_PANC1_no_exclusion.tsv, TableS17_PANC1_grid.tsv,
TableS18_Caki2_grid.tsv and all four with captions in independent_tables.md.
Calls are segments with |log2| >= 0.3; matched = IoU > 0.5 with a reference call of the same direction. HiNT-style supported =
Hi-C calls supported by a reference call with reciprocal 50 % overlap (calls > 2 Mb in parentheses).
"""
import argparse
import os

import pandas as pd

CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX']
W = [3, 5, 7, 9, 11, 21, 41]
DEPMAP = ('depmap_SCC', 'depmap_calls', 'depmap_hits', 'depmap_P', 'depmap_R', 'depmap_bin_SCC')
ARRAY = ('array_SCC', 'arraystate_calls', 'arraystate_hits', 'arraystate_P', 'arraystate_R', 'array_bin_SCC')


def f3(v):
    return f'{float(v):.3f}'


def absent(r):
    """Chromosomes missing from a segmentation: the absent_chromosomes column, else read from the seg file."""
    if 'absent_chromosomes' in r.index and isinstance(r['absent_chromosomes'], str):
        return [] if r['absent_chromosomes'] == 'none' else r['absent_chromosomes'].split(',')
    if not os.path.exists(r['seg']):
        raise SystemExit(f"segmentation {r['seg']} not found; run workflow/20_run_independent.sh first or use --expected")
    present = {l.split('\t')[0] for l in open(r['seg']).readlines()[1:]}
    return [c for c in CHROMS if c not in present]


def pick(d, setting, method, param, lam):
    m = d[(d.setting == setting) & (d.method == method) & (d.param.astype(float) == float(param)) & (d['lambda'].astype(float) == lam)]
    assert len(m) == 1, (setting, method, param, lam)
    return m.iloc[0]


def table7(pe, ck):
    rows = []
    for cell, ref, d, cols, hint in [('PANC-1', 'DepMap WGS', pe, DEPMAP, True), ('PANC-1', 'array-state', pe, ARRAY, False),
                                     ('Caki2', 'DepMap WGS', ck, DEPMAP, True)]:
        for name, method, w in [('baseline', 'none', 0), ('w = 9', 'mean', 9)]:
            r = pick(d, 'complete', method, w, 3.0)
            scc, calls, hits, p, rec, bscc = cols
            rows.append([cell, ref, name + (' †' if absent(r) else ''), f3(r[scc]), str(int(r[calls])), str(int(r[hits])), f3(r[p]), f3(r[rec]),
                         f3(r[bscc]), f"{r['hint_supported']} ({r['hint_large_supported']})" if hint else '—'])
    hdr = ['Cell line', 'Reference', 'Setting', 'Segment SCC', '# calls', '# matched', 'Precision', 'Recall', 'Bin SCC', 'HiNT-style supported (> 2 Mb)']
    cap = ('Table 7. Independent cell lines: baseline (unsmoothed, λ = 3) versus HiNT-RSC (mean filter w = 9, λ = 3), with hg38 problem '
           'regions excluded. Results without problem-region exclusion, including the pre-registered primary comparison (PANC-1 array-state '
           'reference), are given in Supplementary Table S16. † Chromosome 20 is absent from this segmentation.')
    return cap, hdr, rows


def table_s16(pa, pdm):
    rows = []
    for ref, d, pre, hint in [('array-state (pre-registered primary)', pa, 'arraystate', False), ('DepMap WGS', pdm, 'depmap', True)]:
        for name, w in [('baseline', 0), ('w = 9', 9)]:
            m = d[(d.setting == 'complete') & (d.w == w) & (d['lambda'].astype(float) == 3.0)]
            assert len(m) == 1
            r = m.iloc[0]
            rows.append([ref, name, f3(r['SCC']), str(int(r[f'{pre}_calls'])), str(int(r[f'{pre}_hits'])), f3(r[f'{pre}_P']), f3(r[f'{pre}_R']),
                         f3(r['bin_SCC']), r['hint_supported'] if hint else '—'])
    hdr = ['Reference', 'Setting', 'Segment SCC', 'calls', 'matched', 'precision', 'recall', 'bin SCC', 'HiNT-style supported']
    cap = ('Table S16. PANC-1 without hg38 problem-region exclusion: baseline (unsmoothed, λ = 3) versus HiNT-RSC (mean filter w = 9, '
           'λ = 3, both fixed on K562 before the run). The comparison against the array-state reference without exclusion was the '
           'pre-registered primary test; the DepMap WGS reference was added post hoc. Reference calls: DepMap WGS, '
           f"{int(pdm['depmap_ref'].iloc[0])}; array-state, {int(pa['arraystate_ref'].iloc[0])}.")
    return cap, hdr, rows


def grid(d, label, has_array):
    cfg = [('complete', 'none', 0, l, f'unsmoothed, λ = {l:g}') for l in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0)]
    cfg += [('complete', 'mean', w, l, f'mean w = {w}, λ = {l:g}') for l in (3.0, 2.0, 1.0, 0.5) for w in W]
    for mname, lab, g in [('gaussian', 'Gaussian w', W), ('median', 'median w', W), ('savgol', 'Savitzky–Golay w', [5, 7, 9, 11, 21, 41]),
                          ('haar', 'Haar L', [1, 2, 3, 4, 5]), ('sym4', 'sym4 L', [1, 2, 3, 4, 5])]:
        cfg += [('complete', mname, p, 3.0, f'{lab} = {p}, λ = 3') for p in g]
    cfg += [('intra', 'none', 0, 3.0, 'intra-only, unsmoothed, λ = 3')] + [('intra', 'mean', w, 3.0, f'intra-only, mean w = {w}, λ = 3') for w in W]
    rows = []
    for st, m, p, l, lab in cfg:
        r = pick(d, st, m, p, l)
        miss = absent(r)
        row = [lab + (' †' if miss else ''), f3(r['depmap_SCC']), f3(r['depmap_bin_SCC']), str(int(r['depmap_calls'])), str(int(r['depmap_hits'])), f3(r['depmap_P'])]
        if has_array:
            row += [f3(r['array_SCC']), str(int(r['arraystate_hits'])), f3(r['arraystate_P'])]
        rows.append(row + [', '.join(miss) if miss else '—'])
    hdr = ['Setting', 'DepMap SCC', 'DepMap bin SCC', 'calls', 'matched', 'precision'] + (['array SCC', 'array matched', 'array precision'] if has_array else []) + ['absent chromosomes']
    cell = 'PANC-1' if has_array else 'Caki2'
    cap = (f'{label}. Post hoc grids for the independent cell lines ({cell}; hg38 problem regions excluded; complete Hi-C unless stated). '
           'Segment-level SCC, bin-level SCC, CNV calls, matched calls and precision against the DepMap WGS reference'
           + ('; segment-level SCC, matched calls and precision against the ENCODE array-state reference' if has_array else '')
           + '. † and the last column mark segmentations from which BIC-seq2 dropped at least one chromosome.')
    return cap, hdr, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--expected', action='store_true', help='read the frozen results in expected/ instead of output/')
    ap.add_argument('--panc1-excl', help='default: output/panc1/excl/results.tsv')
    ap.add_argument('--caki2-excl', help='default: output/caki2/excl/results.tsv')
    ap.add_argument('--panc1-noexcl-array', help='default: output/panc1/noexcl/results.tsv')
    ap.add_argument('--panc1-noexcl-depmap', help='default: output/panc1/noexcl/results_depmap_PR-Tm9gPE.tsv')
    ap.add_argument('--outdir', default='output/tables')
    args = ap.parse_args()
    defaults = {'panc1_excl': ('output/panc1/excl/results.tsv', 'expected/panc1_excluded_grid.tsv'),
                'caki2_excl': ('output/caki2/excl/results.tsv', 'expected/caki2_excluded_grid.tsv'),
                'panc1_noexcl_array': ('output/panc1/noexcl/results.tsv', 'expected/panc1_noexcl_array.tsv'),
                'panc1_noexcl_depmap': ('output/panc1/noexcl/results_depmap_PR-Tm9gPE.tsv', 'expected/panc1_noexcl_depmap.tsv')}
    for k, (out, exp) in defaults.items():
        if getattr(args, k) is None:
            setattr(args, k, exp if args.expected else out)
    pe, ck, pa, pdm = (pd.read_csv(getattr(args, k), sep='\t') for k in defaults)

    tables = [('Table7_independent_lines.tsv', table7(pe, ck)), ('TableS16_PANC1_no_exclusion.tsv', table_s16(pa, pdm)),
              ('TableS17_PANC1_grid.tsv', grid(pe, 'Table S17', True)), ('TableS18_Caki2_grid.tsv', grid(ck, 'Table S18', False))]
    os.makedirs(args.outdir, exist_ok=True)
    md = []
    for fname, (cap, hdr, rows) in tables:
        pd.DataFrame(rows, columns=hdr).to_csv(os.path.join(args.outdir, fname), sep='\t', index=False)
        md += [f'## {cap}', '', '| ' + ' | '.join(hdr) + ' |', '|' + '---|' * len(hdr)] + ['| ' + ' | '.join(r) + ' |' for r in rows] + ['']
        print('wrote', os.path.join(args.outdir, fname))
    with open(os.path.join(args.outdir, 'independent_tables.md'), 'w') as f:
        f.write('\n'.join(md))
    print('wrote', os.path.join(args.outdir, 'independent_tables.md'))


if __name__ == '__main__':
    main()
