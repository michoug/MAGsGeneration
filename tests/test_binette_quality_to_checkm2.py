"""Tests for scripts/binette_quality_to_checkm2.py."""

import csv
import os
import sys

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"),
)

from binette_quality_to_checkm2 import CHECKM2_HEADER, convert, main

REAL_FIXTURE = "/scratch/gmichoud/test/comparison/binette/336R_mg_25_bins/final_bins_quality_reports.tsv"


def write_binette_tsv(path, rows):
    header = [
        "name",
        "origin",
        "is_original",
        "original_name",
        "completeness",
        "contamination",
        "score",
        "checkm2_model",
        "size",
        "N50",
        "coding_density",
        "contig_count",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class TestConvert:
    def test_maps_known_columns_and_fills_na_for_the_rest(self):
        rows = [
            {
                "name": "binette_bin1",
                "origin": "binette",
                "is_original": "False",
                "original_name": "binette_bin1",
                "completeness": "96.74",
                "contamination": "1.33",
                "score": "94.08",
                "checkm2_model": "Gradient Boost (General Model)",
                "size": "4610981",
                "N50": "4711",
                "coding_density": "0.9",
                "contig_count": "1152",
            }
        ]
        out = convert(rows)
        assert len(out) == 1
        row = out[0]
        assert row["Name"] == "binette_bin1"
        assert row["Completeness"] == "96.74"
        assert row["Contamination"] == "1.33"
        assert row["Completeness_Model_Used"] == "Gradient Boost (General Model)"
        assert row["Coding_Density"] == "0.9"
        assert row["Contig_N50"] == "4711"
        assert row["Genome_Size"] == "4610981"
        assert row["Total_Contigs"] == "1152"
        # Not computed by Binette - must be filled, never fabricated.
        assert row["Translation_Table_Used"] == "NA"
        assert row["Average_Gene_Length"] == "NA"
        assert row["GC_Content"] == "NA"
        assert row["Total_Coding_Sequences"] == "NA"
        assert row["Max_Contig_Length"] == "NA"
        assert row["Additional_Notes"] == "NA"
        assert set(row.keys()) == set(CHECKM2_HEADER)


class TestMain:
    def test_output_header_matches_checkm2_native_schema(self, tmp_path):
        binette_tsv = tmp_path / "final_bins_quality_reports.tsv"
        write_binette_tsv(
            binette_tsv,
            [
                {
                    "name": "binette_bin1",
                    "origin": "binette",
                    "is_original": "False",
                    "original_name": "binette_bin1",
                    "completeness": "96.74",
                    "contamination": "1.33",
                    "score": "94.08",
                    "checkm2_model": "Gradient Boost (General Model)",
                    "size": "4610981",
                    "N50": "4711",
                    "coding_density": "0.9",
                    "contig_count": "1152",
                }
            ],
        )
        out_tsv = tmp_path / "quality_report.tsv"

        rc = main(["--binette-tsv", str(binette_tsv), "--out", str(out_tsv)])

        assert rc == 0
        with open(out_tsv) as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            rows = list(reader)
        assert header == CHECKM2_HEADER
        assert len(rows) == 1
        assert rows[0][header.index("Name")] == "binette_bin1"
        assert rows[0][header.index("Completeness")] == "96.74"
        assert rows[0][header.index("Contamination")] == "1.33"

    def test_against_real_binette_output_if_available(self, tmp_path):
        if not os.path.exists(REAL_FIXTURE):
            return
        out_tsv = tmp_path / "quality_report.tsv"
        rc = main(["--binette-tsv", REAL_FIXTURE, "--out", str(out_tsv)])
        assert rc == 0
        with open(out_tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
        assert len(rows) > 0
        for row in rows:
            comp = float(row["Completeness"])
            cont = float(row["Contamination"])
            assert 0 <= comp <= 100
            assert cont >= 0
            assert row["Name"]
