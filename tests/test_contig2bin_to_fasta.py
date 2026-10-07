"""Tests for scripts/contig2bin_to_fasta.py."""

import os
import sys

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"),
)

from contig2bin_to_fasta import explode, main


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


def write_contig2bin(path, rows):
    with open(path, "w") as fh:
        for contig, bin_id in rows:
            fh.write(f"{contig}\t{bin_id}\n")


class TestExplode:
    def test_splits_bins_and_reports_no_missing(self, tmp_path):
        assembly = tmp_path / "assembly.fa"
        write_fasta(assembly, [("c1", "ACGT"), ("c2", "TTTT"), ("c3", "GGGG")])
        c2b = tmp_path / "c2b.tsv"
        write_contig2bin(c2b, [("c1", "binA"), ("c2", "binA"), ("c3", "binB")])
        outdir = tmp_path / "out"

        n_written, n_missing = explode(str(c2b), str(assembly), str(outdir), min_size=0)

        assert n_written == 2
        assert n_missing == 0
        assert (outdir / "binA.fa").exists()
        assert (outdir / "binB.fa").exists()

    def test_reports_missing_contigs(self, tmp_path):
        assembly = tmp_path / "assembly.fa"
        write_fasta(assembly, [("c1", "ACGT")])
        c2b = tmp_path / "c2b.tsv"
        # c2 is not in the assembly - contig2bin table out of sync.
        write_contig2bin(c2b, [("c1", "binA"), ("c2", "binA")])
        outdir = tmp_path / "out"

        n_written, n_missing = explode(str(c2b), str(assembly), str(outdir), min_size=0)

        assert n_missing == 1

    def test_empty_contig2bin_writes_nothing(self, tmp_path):
        assembly = tmp_path / "assembly.fa"
        write_fasta(assembly, [("c1", "ACGT")])
        c2b = tmp_path / "c2b.tsv"
        c2b.write_text("")
        outdir = tmp_path / "out"

        n_written, n_missing = explode(str(c2b), str(assembly), str(outdir), min_size=0)

        assert n_written == 0
        assert n_missing == 0


class TestMain:
    def test_exits_zero_on_clean_reconstruction(self, tmp_path):
        assembly = tmp_path / "assembly.fa"
        write_fasta(assembly, [("c1", "ACGT")])
        c2b = tmp_path / "c2b.tsv"
        write_contig2bin(c2b, [("c1", "binA")])
        outdir = tmp_path / "out"

        rc = main(
            [
                "--contig2bin",
                str(c2b),
                "--assembly",
                str(assembly),
                "--outdir",
                str(outdir),
            ]
        )

        assert rc == 0

    def test_fails_loudly_on_missing_contigs(self, tmp_path):
        # A contig2bin table referencing contigs absent from the assembly
        # means the two inputs are out of sync - CheckM2 must never score
        # a silently truncated bin built from this mismatch.
        assembly = tmp_path / "assembly.fa"
        write_fasta(assembly, [("c1", "ACGT")])
        c2b = tmp_path / "c2b.tsv"
        write_contig2bin(c2b, [("c1", "binA"), ("ghost_contig", "binA")])
        outdir = tmp_path / "out"

        rc = main(
            [
                "--contig2bin",
                str(c2b),
                "--assembly",
                str(assembly),
                "--outdir",
                str(outdir),
            ]
        )

        assert rc == 1
