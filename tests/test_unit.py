"""Fast tests on synthetic data (no R, BIC-seq2 or .hic files needed)."""
import gzip

import numpy as np
import pytest

from hint_rsc.cli import main
from hint_rsc.covariates import build_covariates
from hint_rsc.evaluate import calls, evaluate, match_calls, read_reference
from hint_rsc.genome import Genome
from hint_rsc.normalize import exclusion_mask
from hint_rsc.segment import log2_ratio, write_bins, write_calls
from hint_rsc.smooth import METHODS, mean_filter, smooth


@pytest.fixture
def tiny(tmp_path):
    sizes = tmp_path / 'tiny.len'
    sizes.write_text('chrA\t120000\nchrB\t60000\nchrY\t50000\n')
    return Genome(str(sizes), 10000)


def test_builtin_genomes():
    g = Genome('hg38')
    assert g.chroms[0] == 'chr1' and g.chroms[-1] == 'chrX' and 'chrY' not in g.length
    assert g.nbins['chr1'] == 248956422 // 50000 + 1
    assert g.builtin_covariates('HindIII').endswith('hg38_HindIII.tsv.gz')
    assert Genome('hg19').builtin_covariates('HindIII') is None
    assert Genome('hg19').builtin_problem_regions() is None and len(g.builtin_problem_regions()) == 3


def test_sizes_file(tiny):
    assert tiny.chroms == ['chrA', 'chrB'] and tiny.nbins == {'chrA': 13, 'chrB': 7}


def test_mean_filter_truncates_at_ends():
    np.testing.assert_allclose(mean_filter(np.arange(1.0, 6.0), 3), [1.5, 2, 3, 4, 4.5])
    np.testing.assert_allclose(mean_filter(np.arange(1.0, 6.0), 1), np.arange(1.0, 6.0))


@pytest.mark.parametrize('method', METHODS)
def test_smoothers_keep_length_and_constant(method):
    x = np.full(200, 0.7)
    p = 3 if method in ('haar', 'sym4') else 9
    y = smooth(x, method, p)
    assert len(y) == 200 and np.allclose(y, 0.7)
    assert smooth(x, method, 0) is x


def test_min_shift_bins(tmp_path):
    config = write_bins({'chrA': np.array([-1.0, 0.0, 2.0])}, str(tmp_path / 'b'), 10000)
    lines = (tmp_path / 'b' / '_chrA.bin').read_text().splitlines()
    assert lines == ['start\tend\tobs\texp\tvar', '0\t9999\t0.0\t1.0\t-1.0', '20000\t29999\t3.0\t1.0\t2.0']
    assert open(config).read().startswith('chr\tpath\nchrA\t')
    assert log2_ratio([-1.0, 0.0, 1.0]) == [0.0, 0.0, 1.0]


def test_cached_segmentation_removed_when_bins_change(tmp_path):
    d = str(tmp_path / 'b')
    write_bins({'chrA': np.array([-1.0, 1.0])}, d)
    (tmp_path / 'b' / 'bicseq_lambda3.seg').write_text('cached')
    write_bins({'chrA': np.array([-1.0, 1.0])}, d)
    assert (tmp_path / 'b' / 'bicseq_lambda3.seg').exists()
    write_bins({'chrA': np.array([-1.0, 2.0])}, d)
    assert not (tmp_path / 'b' / 'bicseq_lambda3.seg').exists()


def test_covariates(tmp_path, tiny):
    seq_a = ('AAGCTT' + 'A' * 9994) + ('G' * 10000) + 'C' * 100000           # sites at bp 1; GC 0 / 1 / 1 ...
    seq_b = 'AT' * 30000
    fa = tmp_path / 'tiny.fa.gz'
    with gzip.open(fa, 'wt') as f:
        f.write('>chrA\n' + seq_a + '\n>chrB\n' + seq_b + '\n')
    bg = tmp_path / 'map.bedGraph'
    bg.write_text('chrA\t0\t5000\t1\nchrA\t5000\t10000\t0.5\n')
    out = tmp_path / 'cov.tsv'
    build_covariates(str(fa), tiny, 'AAGCTT', str(bg), str(out), log=lambda m: None)
    rows = [l.split('\t') for l in out.read_text().splitlines()[1:]]
    assert len(rows) == 13 + 7
    assert rows[0][:3] == ['chrA', '0', '1.0'] and float(rows[0][3]) == pytest.approx(2 / 10000 * 10 / 10, abs=1e-3)
    assert float(rows[1][3]) == 1.0 and float(rows[0][4]) == 0.75 and float(rows[1][4]) == 0.0
    assert float(rows[13][3]) == 0.0


def test_exclusion_mask(tiny):
    # 1-kb resolution, end inclusive: (0, 1400) covers 2 of the 10 kb of bin 0, (0, 2400) covers 3
    m = exclusion_mask([('chrA', 0, 2400), ('chrA', 20000, 30000), ('chrZ', 0, 10)], tiny, 0.25)
    assert m['chrA'][:4].tolist() == [True, False, True, False]
    assert not exclusion_mask([('chrA', 0, 1400)], tiny, 0.25)['chrA'][0]
    assert exclusion_mask([('chrA', 0, 1400)], tiny, 0.2)['chrA'][0]


def test_matching_and_reference_formats(tmp_path, tiny):
    seg = tmp_path / 'hic.seg'
    seg.write_text('chrom\tstart\tend\tbinNum\tobserved\texpected\tlog2.copyRatio\tpvalue\n'
                   'chrA\t0\t39999\t4\t1\t1\t0.5\t0\nchrA\t40000\t119999\t8\t1\t1\t0\t0\nchrB\t0\t59999\t6\t1\t1\t-0.4\t0\n')
    bed = tmp_path / 'ref.bed'
    bed.write_text('chrA\t0\t50000\t0.6\nchrA\t50000\t120000\t0.0\nchrB\t0\t20000\t-1\n')
    ref = read_reference(str(bed), tiny)
    assert ref[0] == ('chrA', 1, 50001, 0.6)
    n_trg, n_ref, hits = match_calls(calls([('chrA', 0, 39999, 0.5), ('chrB', 0, 59999, -0.4)], tiny), calls(ref, tiny, offset=True))
    assert (n_trg, n_ref, hits) == (2, 2, [50000])                    # chrA IoU 0.8 matches; chrB IoU 1/3 does not
    row = evaluate(str(seg), None, ref, tiny)
    assert row['calls'] == 2 and row['matched'] == 1 and row['precision'] == 0.5 and row['recall'] == 0.5
    n = write_calls(str(seg), tiny, str(tmp_path / 'calls.bed'))
    assert n == 2 and (tmp_path / 'calls.bed').read_text().splitlines()[1] == 'chrA\t0\t40000\t0.5000\tgain\t4'


def test_cli_version(capsys):
    with pytest.raises(SystemExit) as e:
        main(['--version'])
    assert e.value.code == 0 and '1.2.0' in capsys.readouterr().out


def test_cli_requires_one_input(tmp_path):
    with pytest.raises(SystemExit, match='one of --hic or --coverage'):
        main(['run', '--genome', 'hg38', '--enzyme', 'HindIII', '-o', str(tmp_path), '--bicseq', __file__])
