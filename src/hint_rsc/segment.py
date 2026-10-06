"""BIC-seq2 input bins (HiNT's per-chromosome min-shift), BIC-seq2 segmentation, CNV calls and bin-level profiles."""
import math
import os
import shutil
import subprocess


def fmt(x):
    """3.0 -> '3', 0.5 -> '0.5' (used in file names and on the BIC-seq2 command line)."""
    return str(int(x)) if float(x).is_integer() else str(float(x))


def write_bins(resid, bin_dir, resolution=50000):
    """Write BIC-seq2 bin files with HiNT-CNV's per-chromosome min-shift: obs = r - min(r), exp = -min(r).

    Bins with a residual of exactly 0 (no GAM fit) are not written. If the bins differ from those already in bin_dir,
    cached segmentations there are removed. Returns the BIC-seq2 config file.
    """
    os.makedirs(bin_dir, exist_ok=True)
    config = os.path.join(bin_dir, '_config.txt')
    changed = False
    lines = ['chr\tpath\n']
    for chrom, values in resid.items():
        path = os.path.abspath(os.path.join(bin_dir, f'_{chrom}.bin'))
        lines.append(f'{chrom}\t{path}\n')
        rows = ['start\tend\tobs\texp\tvar\n']
        min_value = min(values)
        for i, value in enumerate(values):
            obs = value - min_value
            exp = 0 - min_value
            var = obs - exp
            if not var == 0:
                rows.append(f'{i * resolution}\t{(i + 1) * resolution - 1}\t{obs}\t{exp}\t{var}\n')
        changed |= _write_if_changed(path, ''.join(rows))
    changed |= _write_if_changed(config, ''.join(lines))
    if changed:
        for name in os.listdir(bin_dir):
            if name.endswith('.seg'):
                os.remove(os.path.join(bin_dir, name))
    return config


def _write_if_changed(path, text):
    if os.path.exists(path) and open(path).read() == text:
        return False
    with open(path, 'w') as f:
        f.write(text)
    return True


def find_bicseq(path=None):
    for cand in (path, os.environ.get('BICSEQ'), 'NBICseq-seg.pl'):
        if cand and (os.path.exists(cand) or shutil.which(cand)):
            return cand if os.path.exists(cand) else shutil.which(cand)
    raise SystemExit('BIC-seq2 (NBICseq-seg.pl, v0.7.2) not found: install env/bicseq2.yml and set '
                     'BICSEQ="$(conda run -n bicseq2 which NBICseq-seg.pl)", or pass --bicseq')


def bicseq(config, lam, bicseq_pl, tmp_dir, seg=None):
    """Segment with BIC-seq2 at penalty lam (cached: an existing output is reused)."""
    seg = seg or os.path.join(os.path.dirname(config), f'bicseq_lambda{fmt(lam)}.seg')
    if not os.path.exists(seg):
        os.makedirs(tmp_dir, exist_ok=True)
        p = subprocess.run([bicseq_pl, config, seg + '.part', '--lambda', fmt(lam), '--bootstrap', '--tmp', tmp_dir],
                           capture_output=True, text=True)
        if p.returncode != 0 or not os.path.exists(seg + '.part'):
            raise SystemExit(f'BIC-seq2 failed on {config} (lambda {fmt(lam)}):\n{(p.stderr or p.stdout)[-2000:]}')
        os.replace(seg + '.part', seg)
    return seg


def read_seg(seg_file):
    """BIC-seq2 segments: [(chrom, start, end, nbins, log2)]."""
    out = []
    with open(seg_file) as f:
        f.readline()
        for line in f:
            p = line.split()
            out.append((p[0], int(p[1]), int(p[2]), p[3], float(p[6])))
    return out


def write_calls(seg_file, genome, out, threshold=0.3):
    """CNV calls (|log2 copy ratio| >= threshold) as BED: chrom, start, end (half-open), log2 ratio, gain/loss, bins."""
    n = 0
    with open(out, 'w') as f:
        f.write('#chrom\tstart\tend\tlog2_ratio\ttype\tbins\n')
        for chrom, s, e, nb, l2 in read_seg(seg_file):
            if abs(l2) >= threshold:
                f.write(f'{chrom}\t{max(s, 0)}\t{min(e + 1, genome.length[chrom])}\t{l2:.4f}\t{"gain" if l2 > 0 else "loss"}\t{nb}\n')
                n += 1
    return n


def log2_ratio(values):
    """Per-bin log2(obs / exp) of the min-shifted residuals (0 where obs or exp is not positive), as in the paper's bin-level metric."""
    m = min(values)
    out = []
    for v in values:
        obs, exp = v - m, 0 - m
        out.append(math.log2(obs / exp) if obs > 0 and exp > 0 else 0.0)
    return out


def write_profile(resid, fitted, genome, out):
    """Bin-level profile of the (smoothed) residuals over GAM-fitted bins: chrom, start, end, residual, log2_ratio."""
    res = genome.resolution
    with open(out, 'w') as f:
        f.write('chrom\tstart\tend\tresidual\tlog2_ratio\n')
        for chrom, values in resid.items():
            for i, l2 in enumerate(log2_ratio(values)):
                if fitted[chrom][i]:
                    f.write(f'{chrom}\t{i * res}\t{min((i + 1) * res, genome.length[chrom])}\t{values[i]:.6g}\t{l2:.6g}\n')
    return out
