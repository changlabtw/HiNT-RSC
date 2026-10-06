"""Convert an ENCODE/HAIB genotyping-array CNV bed (hg19) into reference seg files for metric.py.

Rules (Supplementary Note S1):
  1. liftOver hg19 -> hg38 (minMatch 0.95); segments that fail are dropped and counted.
  2. Array-state reference calls: amp -> gain, het.del / homo.del -> loss, normal -> neutral.
     Adjacent calls with the same direction and a gap <= 50 kb are merged.
  3. Correlation reference: each segment's logR projected onto bins (bins outside segments = 0).
Outputs (BIC-seq2-like 7-column seg, 1-based starts, read by metric.py with has_offset=True):
  <prefix>.state_calls.seg  value +1 / -1 for calls  (precision/recall, array-state rule)
  <prefix>.logR.seg         value = logR for all segments (correlation; K562-style |logR|>=0.3 sensitivity)

The shipped PANC-1 references (ENCODE experiment ENCSR000AZP) were built with (run from paper/):
  python scripts/build_array_reference.py --bed downloads/ENCFF392GUL.bed.gz --chain downloads/hg19ToHg38.over.chain.gz
         -o data/panc1/array/ENCFF392GUL
  python scripts/build_array_reference.py --bed downloads/ENCFF155BFQ.bed.gz --chain downloads/hg19ToHg38.over.chain.gz
         -o data/panc1/array/ENCFF155BFQ
"""
import argparse
import os
import subprocess
import tempfile

KEEP = {f'chr{i}' for i in range(1, 23)} | {'chrX'}
DIRECTION = {'amp': 1, 'gain': 1, 'het.del': -1, 'homo.del': -1, 'normal': 0}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bed', required=True, help='array CNV bed(.gz): chrom start end state ... logR(last column)')
    ap.add_argument('--chain', required=True, help='UCSC hg19ToHg38.over.chain.gz')
    ap.add_argument('--liftover', default='liftOver', help='UCSC liftOver executable')
    ap.add_argument('--merge-gap', type=int, default=50000)
    ap.add_argument('-o', '--prefix', required=True)
    args = ap.parse_args()

    tmp = tempfile.mkdtemp()
    src = os.path.join(tmp, 'in.bed')
    opener = 'gzip -dc' if args.bed.endswith('.gz') else 'cat'
    rows = subprocess.run(f'{opener} {args.bed}', shell=True, capture_output=True, text=True, check=True).stdout.splitlines()
    with open(src, 'w') as f:
        for i, line in enumerate(rows):
            p = line.split('\t')
            f.write(f'{p[0]}\t{p[1]}\t{p[2]}\t{i}|{p[3]}|{p[-1]}\n')
    lifted, unmapped = os.path.join(tmp, 'out.bed'), os.path.join(tmp, 'unmapped.bed')
    subprocess.run([args.liftover, '-minMatch=0.95', src, args.chain, lifted, unmapped], check=True, capture_output=True)

    segs = []
    for line in open(lifted):
        c, s, e, name = line.split()[:4]
        if c not in KEEP:
            continue
        _, state, logr = name.split('|')
        segs.append((c, int(s), int(e), state, float(logr)))
    n_unmapped = sum(1 for l in open(unmapped) if not l.startswith('#'))
    order = {c: i for i, c in enumerate([f'chr{i}' for i in range(1, 23)] + ['chrX'])}
    segs.sort(key=lambda x: (order[x[0]], x[1]))

    calls = []
    for c, s, e, state, _ in segs:
        d = DIRECTION[state]
        if d == 0:
            continue
        if calls and calls[-1][0] == c and calls[-1][3] == d and s - calls[-1][2] <= args.merge_gap:
            calls[-1][2] = max(calls[-1][2], e)
        else:
            calls.append([c, s, e, d])

    header = 'chrom\tstart\tend\tbinNum\tobserved\texpected\tlog2.copyRatio\n'
    with open(args.prefix + '.state_calls.seg', 'w') as f:
        f.write(header)
        for c, s, e, d in calls:
            f.write(f'{c}\t{s + 1}\t{e + 1}\tNA\tNA\tNA\t{d}\n')
    with open(args.prefix + '.logR.seg', 'w') as f:
        f.write(header)
        for c, s, e, _, logr in segs:
            f.write(f'{c}\t{s + 1}\t{e + 1}\tNA\tNA\tNA\t{logr}\n')
    print(f'input {len(rows)} | lifted {len(segs)} | unmapped {n_unmapped} | array-state calls after merge {len(calls)}')


if __name__ == '__main__':
    main()
