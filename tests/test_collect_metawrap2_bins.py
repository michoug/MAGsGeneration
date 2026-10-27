"""Tests for scripts/collect_metawrap2_bins.py."""

import csv
import os
import sys

import pytest

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"),
)

from collect_metawrap2_bins import collect

HEADER = ["bin", "completeness", "contamination", "GC", "lineage", "N50", "size"]


def setup_run(tmp_path, rows, gunc=None):
    bins = tmp_path / "metawrap_70_10_bins"
    bins.mkdir()
    header = HEADER + (["gunc_css", "gunc_pass"] if gunc is not None else [])
    with open(tmp_path / "stats.tsv", "w") as fh:
        fh.write("\t".join(header) + "\n")
        for i, name in enumerate(rows):
            values = [name, "80.0", "2.0", "0.5", "checkm2", "5000", "2000000"]
            if gunc is not None:
                values += ["0.1", gunc[i]]
            fh.write("\t".join(values) + "\n")
            (bins / f"{name}.fasta").write_text(f">{name}_contig\nACGT\n")
    return bins, tmp_path / "stats.tsv"


def read_table(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class TestCollect:
    def test_renames_with_sample_prefix_in_stats_order(self, tmp_path):
        bins, stats = setup_run(tmp_path, ["bin_002", "bin_001"])
        out = tmp_path / "final_bins"
        kept, dropped = collect(str(bins), str(stats), str(out), str(tmp_path / "t.tsv"), "S1")
        assert (kept, dropped) == (2, 0)
        assert sorted(os.listdir(out)) == ["S1_bin_001.fa", "S1_bin_002.fa"]
        assert (out / "S1_bin_001.fa").read_text().startswith(">bin_002_contig")

    def test_single_bin_set_names_are_normalised(self, tmp_path):
        bins, stats = setup_run(tmp_path, ["S1_concoct_49"])
        out = tmp_path / "final_bins"
        collect(str(bins), str(stats), str(out), str(tmp_path / "t.tsv"), "S1")
        assert os.listdir(out) == ["S1_bin_001.fa"]

    def test_gunc_failures_are_dropped_but_listed(self, tmp_path):
        bins, stats = setup_run(tmp_path, ["bin_001", "bin_002", "bin_003"], gunc=["yes", "no", ""])
        out = tmp_path / "final_bins"
        table = tmp_path / "t.tsv"
        kept, dropped = collect(str(bins), str(stats), str(out), str(table), "S1")
        assert (kept, dropped) == (2, 1)
        assert sorted(os.listdir(out)) == ["S1_bin_001.fa", "S1_bin_002.fa"]
        rows = read_table(table)
        assert [r["kept"] for r in rows] == ["yes", "no", "yes"]
        assert [r["new_name"] for r in rows] == ["S1_bin_001", "", "S1_bin_002"]

    def test_keep_chimeric(self, tmp_path):
        bins, stats = setup_run(tmp_path, ["bin_001"], gunc=["no"])
        out = tmp_path / "final_bins"
        kept, dropped = collect(
            str(bins), str(stats), str(out), str(tmp_path / "t.tsv"), "S1", drop_chimeric=False
        )
        assert (kept, dropped) == (1, 0)

    def test_header_only_stats_gives_empty_outputs(self, tmp_path):
        bins, stats = setup_run(tmp_path, [])
        out = tmp_path / "final_bins"
        table = tmp_path / "t.tsv"
        assert collect(str(bins), str(stats), str(out), str(table), "S1") == (0, 0)
        assert os.listdir(out) == []
        assert read_table(table) == []

    def test_missing_bin_file_is_an_error(self, tmp_path):
        bins, stats = setup_run(tmp_path, ["bin_001"])
        os.remove(bins / "bin_001.fasta")
        with pytest.raises(FileNotFoundError):
            collect(str(bins), str(stats), str(tmp_path / "o"), str(tmp_path / "t.tsv"), "S1")
