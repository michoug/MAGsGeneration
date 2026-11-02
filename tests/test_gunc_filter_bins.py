"""Tests for scripts/gunc_filter_bins.py."""

import csv
import os
import sys

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"),
)

from gunc_filter_bins import decide, filter_bins, main


def q(cont, comp=90):
    return {"Name": "b", "Completeness": str(comp), "Contamination": str(cont)}


class TestDecide:
    def test_pass(self):
        assert decide(q(5), {"pass.GUNC": "True"}, 10) == (True, "pass")

    def test_chimeric(self):
        assert decide(q(5), {"pass.GUNC": "False"}, 10) == (False, "gunc_chimeric")

    def test_minor_chimeric_kept(self):
        g = {"pass.GUNC": "False", "contamination_portion": "0.05"}
        assert decide(q(5), g, 10) == (True, "pass;gunc_minor_chimeric")
        g = {"pass.GUNC": "False", "contamination_portion": "0.06"}
        assert decide(q(5), g, 10) == (False, "gunc_chimeric")

    def test_minor_chimeric_still_needs_low_contamination(self):
        g = {"pass.GUNC": "False", "contamination_portion": "0.02"}
        assert decide(q(12), g, 10) == (False, "contamination>10")

    def test_portion_cutoff_zero_is_strict(self):
        g = {"pass.GUNC": "False", "contamination_portion": "0.02"}
        assert decide(q(5), g, 10, max_chimeric_portion=0) == (False, "gunc_chimeric")

    def test_contamination_ceiling_inclusive(self):
        assert decide(q(10), {"pass.GUNC": "True"}, 10)[0] is True
        assert decide(q(10.01), {"pass.GUNC": "True"}, 10) == (False, "contamination>10")

    def test_both_reasons(self):
        assert decide(q(15), {"pass.GUNC": "False"}, 10) == (
            False,
            "gunc_chimeric;contamination>10",
        )

    def test_gunc_undecided_or_missing_is_kept(self):
        assert decide(q(5), {"pass.GUNC": "nan"}, 10)[0] is True
        assert decide(q(5), None, 10)[0] is True


def setup(tmp_path):
    bins = tmp_path / "final_bins"
    bins.mkdir()
    for name in ["binette_bin1", "binette_bin2", "binette_bin3", "binette_bin4"]:
        (bins / f"{name}.fa").write_text(f">{name}_c1\nACGT\n")
    (bins / ".snakemake_timestamp").write_text("")
    quality = tmp_path / "quality_report.tsv"
    quality.write_text(
        "Name\tCompleteness\tContamination\n"
        "binette_bin1\t95\t2\n"
        "binette_bin2\t80\t15\n"
        "binette_bin3\t90\t3\n"
        "binette_bin4\t70\t8\n"
    )
    gunc = tmp_path / "GUNC.maxCSS_level.tsv"
    gunc.write_text(
        "genome\tclade_separation_score\tcontamination_portion\tpass.GUNC\n"
        "binette_bin1\t0.1\t0.05\tTrue\n"
        "binette_bin2\t0.2\t0.1\tTrue\n"
        "binette_bin3\t0.9\t0.4\tFalse\n"
        "binette_bin4\tnan\tnan\tnan\n"
    )
    return bins, quality, gunc


class TestFilterBins:
    def test_copies_kept_bins_and_writes_table(self, tmp_path):
        bins, quality, gunc = setup(tmp_path)
        out = tmp_path / "gunc_bins"
        table = tmp_path / "gunc_filter.tsv"
        kept, total = filter_bins(str(bins), str(quality), str(gunc), str(out), str(table))
        assert (kept, total) == (2, 4)
        assert sorted(os.listdir(out)) == ["binette_bin1.fa", "binette_bin4.fa"]
        rows = {r["bin"]: r for r in csv.DictReader(open(table), delimiter="\t")}
        assert rows["binette_bin2"]["reason"] == "contamination>10"
        assert rows["binette_bin3"]["reason"] == "gunc_chimeric"
        assert rows["binette_bin4"]["kept"] == "yes"
        assert rows["binette_bin3"]["gunc_css"] == "0.9"

    def test_main_max_contamination(self, tmp_path):
        bins, quality, gunc = setup(tmp_path)
        out = tmp_path / "gunc_bins"
        table = tmp_path / "t.tsv"
        args = ["--bins", str(bins), "--quality", str(quality), "--gunc", str(gunc),
                "--outdir", str(out), "--table", str(table), "--max-contamination", "20"]
        assert main(args) == 0
        assert sorted(os.listdir(out)) == ["binette_bin1.fa", "binette_bin2.fa", "binette_bin4.fa"]

    def test_rerun_removes_stale_bins(self, tmp_path):
        bins, quality, gunc = setup(tmp_path)
        out = tmp_path / "gunc_bins"
        table = tmp_path / "t.tsv"
        base = ["--bins", str(bins), "--quality", str(quality), "--gunc", str(gunc),
                "--outdir", str(out), "--table", str(table)]
        main(base + ["--max-contamination", "20"])
        assert "binette_bin2.fa" in os.listdir(out)
        main(base + ["--max-contamination", "10"])
        assert sorted(os.listdir(out)) == ["binette_bin1.fa", "binette_bin4.fa"]

    def test_empty_bins(self, tmp_path):
        bins = tmp_path / "empty"
        bins.mkdir()
        quality = tmp_path / "q.tsv"
        quality.write_text("Name\tCompleteness\tContamination\n")
        gunc = tmp_path / "g.tsv"
        gunc.write_text("genome\tpass.GUNC\n")
        out = tmp_path / "out"
        assert filter_bins(str(bins), str(quality), str(gunc), str(out), str(tmp_path / "t.tsv")) == (0, 0)
        assert os.listdir(out) == []
