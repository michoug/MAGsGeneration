"""Tests for scripts/read_counts_mqc.py."""

import json
import os
import sys

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"),
)

from read_counts_mqc import counts, main, sample_name

TRIM = {
    "pair_validation": {
        "pairs_analyzed": 1000,
        "pairs_removed": 10,
        "pairs_removed_n": 5,
        "pairs_removed_too_long": 0,
    }
}
FASTP = {
    "summary": {
        "before_filtering": {"total_reads": 1960},
        "after_filtering": {"total_reads": 1900},
    }
}


class TestCounts:
    def test_pairs_per_step(self):
        c = counts(TRIM, FASTP)
        assert c["raw_pairs"] == 1000
        assert c["trimmed_pairs"] == 985
        assert c["human_removed_pairs"] == 985 - 980
        assert c["fastp_removed_pairs"] == 980 - 950
        assert c["final_pairs"] == 950
        assert c["pct_retained"] == 95.0

    def test_pct_human(self):
        assert counts(TRIM, FASTP)["pct_human"] == round(100 * 5 / 985, 3)

    def test_missing_optional_counters(self):
        trim = {"pair_validation": {"pairs_analyzed": 100, "pairs_removed": 0}}
        fastp = {
            "summary": {
                "before_filtering": {"total_reads": 200},
                "after_filtering": {"total_reads": 200},
            }
        }
        c = counts(trim, fastp)
        assert c["trimmed_pairs"] == 100
        assert c["human_removed_pairs"] == 0

    def test_empty_sample(self):
        trim = {"pair_validation": {"pairs_analyzed": 0, "pairs_removed": 0}}
        fastp = {
            "summary": {
                "before_filtering": {"total_reads": 0},
                "after_filtering": {"total_reads": 0},
            }
        }
        c = counts(trim, fastp)
        assert c["pct_human"] == 0
        assert c["pct_retained"] == 0


class TestMain:
    def test_sample_name(self):
        assert sample_name("/x/TrimFinal/Fastp/S1_L001_fastp.json") == "S1_L001"

    def test_writes_mqc_table(self, tmp_path):
        trim_dir = tmp_path / "Trim"
        trim_dir.mkdir()
        (trim_dir / "S1_R2_001.fastq.gz_trimming_report.json").write_text(json.dumps(TRIM))
        fastp = tmp_path / "S1_fastp.json"
        fastp.write_text(json.dumps(FASTP))
        out = tmp_path / "read_counts_mqc.tsv"

        assert main(["--fastp", str(fastp), "--trim-dir", str(trim_dir), "--out", str(out)]) == 0
        lines = out.read_text().splitlines()
        assert lines[0].startswith("# id: 'read_counts'")
        data = [l for l in lines if not l.startswith("#")]
        assert data[0].split("\t")[0] == "Sample"
        assert data[1].split("\t")[:2] == ["S1", "1000"]
        assert data[1].split("\t")[-2] == "950"
