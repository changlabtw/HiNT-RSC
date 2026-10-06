"""Re-score cached run_validation.py segmentations against the DepMap WGS reference, with repo metrics and HiNT-paper metrics.

DepMap reference rule (Supplementary Note S1): log2(SegmentMean) (floor -5), calls |log2| >= 0.3.
HiNT-paper "supported" rule (HiNT Fig. S5C): a Hi-C call is supported if a reference call of the same direction
overlaps > 50% of the Hi-C segment and > 50% of the reference segment.
Writes <outdir>/depmap_<profile>.seg and <outdir>/results_depmap_<profile>.tsv (PANC-1 without problem-region exclusion:
DepMap rows of Supplementary Table S16). Genome: hg38 (data/ref_genome/hg38).
"""
import argparse
import csv
import glob
import math
import os

if os.path.basename(os.path.normpath(os.environ.get('HINT_RSC_REF', ''))) != 'hg38':
    os.environ['HINT_RSC_REF'] = 'data/ref_genome/hg38'  # PANC-1 and Caki2 are analysed on hg38

import pandas as pd  # noqa: E402

from load import chr_length, load_seg_pair  # noqa: E402
from run_validation import score  # noqa: E402


def depmap_to_seg(csv_file, profile, out_seg):
    """Write one DepMap profile (OmicsCNSegmentsProfile.csv rows) as a 7-column seg file with log2 ratios."""
    with open(out_seg, 'w') as f:
        f.write('chrom\tstart\tend\tbinNum\tobserved\texpected\tlog2.copyRatio\n')
        for r in csv.DictReader(open(csv_file)):
            if r['ProfileID'] != profile:
                continue
            chrom = 'chr' + r['Chromosome']
            if chrom not in chr_length:
                continue
            m = float(r['SegmentMean'])
            log2 = math.log2(m) if m > 0 else -5.0
            f.write(f"{chrom}\t{int(r['Start']) + 1}\t{int(r['End']) + 1}\tNA\tNA\tNA\t{max(log2, -5.0)}\n")
    return out_seg


def hint_supported(seg_file, ref_seg, large=2_000_000):
    trg = load_seg_pair(seg_file, threshold=0.3)
    ref = load_seg_pair(ref_seg, threshold=0.3, has_offset=True)
    n = sup = n_large = sup_large = 0
    for chrom, calls in trg.items():
        for s, e, d in calls:
            ok = False
            for rs, re_, rd in ref[chrom]:
                ov = max(0, min(e, re_) - max(s, rs))
                if rd == d and ov > 0.5 * (e - s) and ov > 0.5 * (re_ - rs):
                    ok = True
                    break
            n += 1
            sup += ok
            if e - s > large:
                n_large += 1
                sup_large += ok
    return sup, n, sup_large, n_large


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--depmap-csv', default='data/panc1/depmap/PANC1_PR-Tm9gPE_CNsegments.csv',
                    help='DepMap OmicsCNSegmentsProfile.csv or a per-profile extract (default: %(default)s)')
    ap.add_argument('--profile', default='PR-Tm9gPE', help='DepMap ProfileID (default: %(default)s, PANC-1 WGS)')
    ap.add_argument('--outdir', default='output/panc1/noexcl',
                    help='run_validation.py output directory with cached segmentations (default: %(default)s)')
    args = ap.parse_args()

    ref_seg = depmap_to_seg(args.depmap_csv, args.profile, os.path.join(args.outdir, f'depmap_{args.profile}.seg'))
    refs = {'corr': ref_seg, 'calls': {'depmap': ref_seg}}
    rows = []
    for setting in ('complete', 'intra'):
        for seg_file in sorted(glob.glob(os.path.join(args.outdir, setting, 'w*', 'w*_lambda*.seg'))):
            w = int(os.path.basename(os.path.dirname(seg_file))[1:])
            lam = os.path.basename(seg_file).split('lambda')[1][:-4]
            row = {'setting': setting, 'w': w, 'lambda': lam, **score(seg_file, os.path.dirname(seg_file), refs)}
            sup, n, sup_l, n_l = hint_supported(seg_file, ref_seg)
            row.update({'hint_supported': f'{sup}/{n}', 'hint_supported_frac': round(sup / n, 3) if n else None,
                        'hint_large_supported': f'{sup_l}/{n_l}', 'hint_large_frac': round(sup_l / n_l, 3) if n_l else None})
            rows.append(row)
    df = pd.DataFrame(rows).sort_values(['setting', 'w', 'lambda'])
    out = os.path.join(args.outdir, f'results_depmap_{args.profile}.tsv')
    df.to_csv(out, sep='\t', index=False)
    print(df.to_string(index=False))
    print('wrote', out)


if __name__ == '__main__':
    main()
