"""Per-sample bin counts by quality tier, as a MultiQC custom-content table.

Reads the per-sample CheckM2 quality reports of the refined bins
(<sample>/checkm2/quality_report.tsv) and counts the bins in each
completeness / contamination tier. Optionally also writes the binner
selection report (binner_selection_report.tsv) as a second table.
"""

import argparse
import csv
import os
import sys

TIERS = [
    ("mq50", 50, 10),
    ("mq70", 70, 10),
    ("hq90", 90, 5),
    ("nc95", 95, 5),
]
MAG_HEADER = """# id: 'mag_counts'
# section_name: 'Refined bins per sample'
# description: 'Refined bins per sample by CheckM2 tier: MQ >=50% / <10%, MQ >=70% / <10%, HQ >=90% / <5%, NC >=95% / <5% (completeness / contamination).'
# plot_type: 'table'
# pconfig:
#     id: 'mag_counts_table'
#     title: 'Refined bins per sample'
"""
BINNER_HEADER = """# id: 'binner_selection'
# section_name: 'Binner selection'
# description: 'CheckM2-based ranking of the binners (lower score is better); selected binners feed the refinement.'
# plot_type: 'table'
# pconfig:
#     id: 'binner_selection_table'
#     title: 'Binner selection'
"""


def sample_of(report_path):
    """<sample> from .../<sample>/checkm2/quality_report.tsv."""
    return os.path.basename(os.path.dirname(os.path.dirname(report_path)))


def tier_counts(rows):
    out = {"bins": len(rows)}
    for name, comp, cont in TIERS:
        out[name] = sum(
            float(r["Completeness"]) >= comp and float(r["Contamination"]) < cont
            for r in rows
        )
    return out


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_mag_table(reports, out):
    columns = ["bins"] + [name for name, _, _ in TIERS]
    out.write(MAG_HEADER)
    out.write("\t".join(["Sample"] + columns) + "\n")
    for path in sorted(reports, key=sample_of):
        counts = tier_counts(read_tsv(path))
        out.write("\t".join([sample_of(path)] + [str(counts[c]) for c in columns]) + "\n")


def write_binner_table(report, out):
    rows = read_tsv(report)
    columns = [c for c in rows[0] if c != "binner"] if rows else []
    out.write(BINNER_HEADER)
    out.write("\t".join(["Binner"] + columns) + "\n")
    for r in rows:
        out.write("\t".join([r["binner"]] + [r[c] for c in columns]) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--binner-report", help="binner_selection_report.tsv")
    parser.add_argument("--binner-out", help="output for the binner table")
    args = parser.parse_args(argv)

    with open(args.out, "w", newline="") as out:
        write_mag_table(args.reports, out)
    if args.binner_report:
        with open(args.binner_out, "w", newline="") as out:
            write_binner_table(args.binner_report, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
