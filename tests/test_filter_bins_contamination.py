"""Tests for scripts/filter_bins_contamination.py."""

import csv
import os
import sys

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"),
)

from filter_bins_contamination import filter_bins, keep, main

REPORT = (
    "Name\tCompleteness\tContamination\tGenome_Size\n"
    "S1_bin1\t95\t2.5\t3000000\n"
    "S1_bin2\t80\t10\t2000000\n"
    "S2_bin1\t88\t10.01\t2500000\n"
    "S2_bin2\t72\t19.5\t4000000\n"
)


def setup(tmp_path):
    bins = tmp_path / "renamed_mags"
    bins.mkdir()
    for name in ["S1_bin1", "S1_bin2", "S2_bin1", "S2_bin2"]:
        (bins / f"{name}.fa").write_text(f">{name}_c1\nACGT\n")
    report = tmp_path / "quality_report.tsv"
    report.write_text(REPORT)
    return bins, report


class TestKeep:
    def test_ceiling_is_inclusive(self):
        assert keep({"Contamination": "10"}, 10)
        assert not keep({"Contamination": "10.01"}, 10)


class TestFilterBins:
    def test_copies_and_filters_report(self, tmp_path):
        bins, report = setup(tmp_path)
        out = tmp_path / "filtered_mags"
        out_report = tmp_path / "filtered.tsv"
        assert filter_bins(str(bins), str(report), str(out), str(out_report)) == (2, 4)
        assert sorted(os.listdir(out)) == ["S1_bin1.fa", "S1_bin2.fa"]
        rows = list(csv.DictReader(open(out_report), delimiter="\t"))
        assert [r["Name"] for r in rows] == ["S1_bin1", "S1_bin2"]
        assert list(rows[0]) == ["Name", "Completeness", "Contamination", "Genome_Size"]

    def test_main_custom_ceiling(self, tmp_path):
        bins, report = setup(tmp_path)
        out = tmp_path / "out"
        out_report = tmp_path / "r.tsv"
        args = ["--bins", str(bins), "--report", str(report), "--outdir", str(out),
                "--out-report", str(out_report), "--max-contamination", "20"]
        assert main(args) == 0
        assert len(os.listdir(out)) == 4

    def test_empty_report(self, tmp_path):
        bins = tmp_path / "b"
        bins.mkdir()
        report = tmp_path / "q.tsv"
        report.write_text("Name\tCompleteness\tContamination\n")
        out_report = tmp_path / "r.tsv"
        assert filter_bins(str(bins), str(report), str(tmp_path / "o"), str(out_report)) == (0, 0)
        assert out_report.read_text() == "Name\tCompleteness\tContamination\n"
