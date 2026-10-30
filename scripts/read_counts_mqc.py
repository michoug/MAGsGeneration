"""Read pairs kept at each read-processing step, as a MultiQC custom-content table.

Takes the R2 Trim Galore (v2+) JSON report of each sample, whose
`pair_validation` block gives the pairs before and after trimming, and the
fastp JSON report, whose `before_filtering` / `after_filtering` read counts
give the pairs after human removal and the final pairs. Pairs removed as
human are the trimmed pairs minus the pairs fastp received.
"""

import argparse
import json
import os
import sys

HEADER = """# id: 'read_counts'
# section_name: 'Read pairs per step'
# description: 'Read pairs after each processing step: Trim Galore, human read removal (strobealign vs hg38), fastp.'
# plot_type: 'table'
# pconfig:
#     id: 'read_counts_table'
#     title: 'Read pairs per step'
"""
COLUMNS = [
    "Sample",
    "raw_pairs",
    "trimmed_pairs",
    "human_removed_pairs",
    "pct_human",
    "fastp_removed_pairs",
    "final_pairs",
    "pct_retained",
]


def counts(trim_report, fastp_report):
    """Pair counts for one sample from its R2 Trim Galore and fastp JSON reports."""
    pv = trim_report["pair_validation"]
    raw = pv["pairs_analyzed"]
    trimmed = raw - pv["pairs_removed"] - pv.get("pairs_removed_n", 0)
    trimmed -= pv.get("pairs_removed_too_long", 0)
    before = fastp_report["summary"]["before_filtering"]["total_reads"] // 2
    final = fastp_report["summary"]["after_filtering"]["total_reads"] // 2
    human = trimmed - before
    return {
        "raw_pairs": raw,
        "trimmed_pairs": trimmed,
        "human_removed_pairs": human,
        "pct_human": round(100 * human / trimmed, 3) if trimmed else 0,
        "fastp_removed_pairs": before - final,
        "final_pairs": final,
        "pct_retained": round(100 * final / raw, 2) if raw else 0,
    }


def sample_name(fastp_json):
    return os.path.basename(fastp_json).removesuffix("_fastp.json")


def write_table(rows, out):
    out.write(HEADER)
    out.write("\t".join(COLUMNS) + "\n")
    for sample, row in rows:
        out.write("\t".join([sample] + [str(row[c]) for c in COLUMNS[1:]]) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastp", nargs="+", required=True, help="{sample}_fastp.json")
    parser.add_argument("--trim-dir", required=True, help="Trim Galore output folder")
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    rows = []
    for fastp_json in sorted(args.fastp):
        sample = sample_name(fastp_json)
        trim_json = os.path.join(
            args.trim_dir, f"{sample}_R2_001.fastq.gz_trimming_report.json"
        )
        with open(trim_json) as fh:
            trim_report = json.load(fh)
        with open(fastp_json) as fh:
            fastp_report = json.load(fh)
        rows.append((sample, counts(trim_report, fastp_report)))
    with open(args.out, "w") as out:
        write_table(rows, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
