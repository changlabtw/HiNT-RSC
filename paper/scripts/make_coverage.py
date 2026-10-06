"""Compute 1D Hi-C coverage profiles (complete = intra + inter, and intra-only) from a .hic file with the load.py functions.

Output columns: chrom, bin_num, rowsum_complete, rowsum_intra (one row per 50-kb bin).
The shipped PANC-1 and Caki2 tables were built with (run from paper/):
  python scripts/make_coverage.py --hic downloads/ENCFF545SGZ.hic -o data/panc1/PANC1_coverage_50kb.txt
  python scripts/make_coverage.py --hic downloads/ENCFF356GEH.hic -o data/caki2/Caki2_coverage_50kb.txt
"""
import argparse
import os
import time


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--hic', required=True)
    ap.add_argument('--ref', default='data/ref_genome/hg38',
                    help='genome folder containing <name>.len (sets HINT_RSC_REF; default: %(default)s)')
    ap.add_argument('--resolution', type=int, default=50000)
    ap.add_argument('-o', '--output', required=True)
    args = ap.parse_args()

    os.environ['HINT_RSC_REF'] = args.ref
    from load import chr_length, generate_coverage_profile_intra, generate_coverage_profile_intra_inter

    with open(args.output, 'w') as f:
        f.write('chrom\tbin_num\trowsum_complete\trowsum_intra\n')
        for chrom in chr_length:
            t = time.time()
            complete = generate_coverage_profile_intra_inter(args.hic, trg_chr=chrom, resolution=args.resolution)
            intra = generate_coverage_profile_intra(args.hic, trg_chr=chrom, resolution=args.resolution)
            for i in range(len(complete)):
                f.write(f'{chrom}\t{i}\t{complete[i]}\t{intra[i]}\n')
            print(f'{chrom}: {len(complete)} bins, {time.time() - t:.0f}s', flush=True)


if __name__ == '__main__':
    main()
