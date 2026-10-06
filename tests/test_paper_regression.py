"""hint-rsc reproduces the paper (paper/expected/) on the shipped PANC-1 and Caki2 inputs.

Needs Rscript with mgcv and BIC-seq2 (NBICseq-seg.pl on PATH or $BICSEQ); skipped otherwise. About 3 minutes.
The coverage test also needs the ENCODE PANC-1 .hic file: HINT_RSC_TEST_HIC=/path/to/ENCFF545SGZ.hic.
"""
import csv
import gzip
import math
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from hint_rsc.cli import main
from hint_rsc.genome import resource

PAPER = Path(__file__).resolve().parents[1] / 'paper'


def _tools():
    r = shutil.which('Rscript')
    ok_r = r and subprocess.run([r, '-e', 'library(mgcv); library(optparse)'], capture_output=True).returncode == 0
    bic = os.environ.get('BICSEQ') or shutil.which('NBICseq-seg.pl')
    return ok_r and bic


needs_tools = pytest.mark.skipif(not _tools(), reason='needs Rscript (mgcv, optparse) and BIC-seq2')

LINES = {'panc1': ('data/panc1/PANC1_coverage_50kb.txt', 'data/panc1/depmap/PANC1_PR-Tm9gPE_CNsegments.csv', 'PR-Tm9gPE',
                   'expected/panc1_excluded_grid.tsv'),
         'caki2': ('data/caki2/Caki2_coverage_50kb.txt', 'data/caki2/depmap/Caki2_PR-g3gtC7_CNsegments.csv', 'PR-g3gtC7',
                   'expected/caki2_excluded_grid.tsv')}
METRICS = {'seg_PCC': 'depmap_PCC', 'seg_SCC': 'depmap_SCC', 'bin_PCC': 'depmap_bin_PCC', 'bin_SCC': 'depmap_bin_SCC',
           'calls': 'depmap_calls', 'ref_calls': 'depmap_ref', 'matched': 'depmap_hits', 'precision': 'depmap_P',
           'recall': 'depmap_R', 'matched_by_ref_size': 'depmap_size'}


def depmap_bed(csv_file, profile, out):
    """DepMap segments as BED with log2(SegmentMean) (floor -5), as in the paper."""
    with open(out, 'w') as f:
        for r in csv.DictReader(open(csv_file)):
            if r['ProfileID'] == profile:
                m = float(r['SegmentMean'])
                f.write(f"chr{r['Chromosome']}\t{r['Start']}\t{r['End']}\t{max(math.log2(m) if m > 0 else -5.0, -5.0)}\n")
    return str(out)


def expected(line, method, param, lam):
    rows = csv.DictReader(open(PAPER / LINES[line][3]), delimiter='\t')
    return next(r for r in rows if r['setting'] == 'complete' and r['method'] == method and int(r['param']) == param
                and float(r['lambda']) == lam)


def same(got, exp):
    for k, col in METRICS.items():
        g, e = got[k], exp[col]
        assert (str(g) == e) if k == 'matched_by_ref_size' else float(g) == float(e), f'{k}: {g} != {e}'


def read_tsv(path):
    return list(csv.DictReader(open(path), delimiter='\t'))


def test_bundled_covariates_match_paper():
    for name, paper in [('hg19_MboI.tsv.gz', 'hg19_MboI_cov.txt'), ('hg38_HindIII.tsv.gz', 'hg38_HindIII_cov.txt')]:
        assert gzip.open(resource('covariates', name), 'rt').read() == (PAPER / 'data/covariates' / paper).read_text()


@needs_tools
@pytest.mark.parametrize('line', ['panc1', 'caki2'])
def test_sweep_reproduces_grid(tmp_path, line):
    cov, dm, profile, _ = LINES[line]
    ref = depmap_bed(PAPER / dm, profile, tmp_path / 'depmap.bed')
    main(['sweep', '--coverage', str(PAPER / cov), '--genome', 'hg38', '--enzyme', 'HindIII', '--windows', '0,9',
          '--lambdas', '1,3', '--workers', '4', '--reference', ref, '-o', str(tmp_path / 'sweep')])
    for row in read_tsv(tmp_path / 'sweep' / 'sweep.tsv'):
        w = int(row['window'])
        same(row, expected(line, 'none' if w == 0 else 'mean', w, float(row['lambda'])))


@needs_tools
@pytest.mark.parametrize('method,param', [('gaussian', 9), ('median', 9), ('savgol', 9), ('haar', 3), ('sym4', 3)])
def test_filters_reproduce_grid(tmp_path_factory, method, param):
    out = tmp_path_factory.getbasetemp() / 'panc1_filters'          # shared: the GAM is fitted once
    ref = depmap_bed(PAPER / LINES['panc1'][1], LINES['panc1'][2], tmp_path_factory.getbasetemp() / 'panc1_depmap.bed')
    main(['run', '--coverage', str(PAPER / LINES['panc1'][0]), '--genome', 'hg38', '--enzyme', 'HindIII', '--smoother', method,
          '-w', str(param), '--reference', ref, '-o', str(out)])
    same(read_tsv(out / 'metrics.tsv')[0], expected('panc1', method, param, 3))


@pytest.mark.skipif(not os.environ.get('HINT_RSC_TEST_HIC'), reason='set HINT_RSC_TEST_HIC to the PANC-1 .hic (ENCFF545SGZ)')
def test_coverage_from_hic(tmp_path):
    out = tmp_path / 'cov.tsv'
    main(['coverage', '--hic', os.environ['HINT_RSC_TEST_HIC'], '--genome', 'hg38', '--chroms', 'chr21', '-o', str(out)])
    exp = [l for l in (PAPER / LINES['panc1'][0]).read_text().splitlines() if l.startswith('chr21\t')]
    assert out.read_text().splitlines()[1:] == exp
