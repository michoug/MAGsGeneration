"""Tests for scripts/mag_counts_mqc.py."""

import os
import sys

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"),
)

from mag_counts_mqc import main, sample_of, tier_counts

HEADER = "Name\tCompleteness\tContamination\n"


def bins(*pairs):
    return [{"Completeness": str(c), "Contamination": str(x)} for c, x in pairs]


def write_report(root, sample, pairs):
    d = root / sample / "checkm2"
    d.mkdir(parents=True)
    path = d / "quality_report.tsv"
    path.write_text(HEADER + "".join(f"b{i}\t{c}\t{x}\n" for i, (c, x) in enumerate(pairs)))
    return path


class TestTierCounts:
    def test_tiers(self):
        c = tier_counts(bins((96, 1), (91, 4), (75, 8), (55, 9), (40, 1), (99, 12)))
        assert c == {"bins": 6, "mq50": 4, "mq70": 3, "hq90": 2, "nc95": 1}

    def test_thresholds_inclusive_completeness_exclusive_contamination(self):
        c = tier_counts(bins((90, 4.99), (90, 5), (50, 9.99), (50, 10)))
        assert c["hq90"] == 1
        assert c["mq50"] == 3

    def test_no_bins(self):
        assert tier_counts([]) == {"bins": 0, "mq50": 0, "mq70": 0, "hq90": 0, "nc95": 0}


class TestMain:
    def test_sample_of(self):
        assert sample_of("/r/Bins/intermediate/S1/checkm2/quality_report.tsv") == "S1"

    def test_mag_table(self, tmp_path):
        r1 = write_report(tmp_path, "S2", [(95, 1), (60, 2)])
        r2 = write_report(tmp_path, "S1", [(80, 3)])
        out = tmp_path / "mag_counts_mqc.tsv"
        assert main(["--reports", str(r1), str(r2), "--out", str(out)]) == 0
        data = [l.split("\t") for l in out.read_text().splitlines() if not l.startswith("#")]
        assert data[0] == ["Sample", "bins", "mq50", "mq70", "hq90", "nc95"]
        assert data[1] == ["S1", "1", "1", "1", "0", "0"]
        assert data[2] == ["S2", "2", "2", "1", "1", "1"]

    def test_binner_table(self, tmp_path):
        r = write_report(tmp_path, "S1", [(95, 1)])
        rep = tmp_path / "binner_selection_report.tsv"
        rep.write_text("binner\tscore\tselected\ncomebin\t1.3\tTrue\nrosella\t6.0\tFalse\n")
        out = tmp_path / "mags.tsv"
        bout = tmp_path / "binner_selection_mqc.tsv"
        main(["--reports", str(r), "--out", str(out), "--binner-report", str(rep), "--binner-out", str(bout)])
        lines = bout.read_text().splitlines()
        assert lines[0] == "# id: 'binner_selection'"
        data = [l.split("\t") for l in lines if not l.startswith("#")]
        assert data[0] == ["Binner", "score", "selected"]
        assert data[1] == ["comebin", "1.3", "True"]
