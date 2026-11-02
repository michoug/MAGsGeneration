"""Keep the bins of a CheckM2 report whose contamination is at most a ceiling.

Copies <bins>/<Name><suffix> for every row of the CheckM2 quality report with
Contamination <= --max-contamination to --outdir, and writes the same report
restricted to those rows to --out-report.
"""

import argparse
import csv
import os
import shutil
import sys


def keep(row, max_contamination):
    return float(row["Contamination"]) <= max_contamination


def filter_bins(bins_dir, report, outdir, out_report, max_contamination=10.0, suffix=".fa"):
    """Copy kept bins and write the filtered report; return (kept, total)."""
    with open(report, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or ["Name", "Completeness", "Contamination"]
        rows = list(reader)
    kept = [r for r in rows if keep(r, max_contamination)]
    os.makedirs(outdir, exist_ok=True)
    for r in kept:
        fn = r["Name"] + suffix
        shutil.copy(os.path.join(bins_dir, fn), os.path.join(outdir, fn))
    with open(out_report, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(kept)
    return len(kept), len(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bins", required=True)
    parser.add_argument("--report", required=True, help="CheckM2 quality_report.tsv")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--out-report", required=True)
    parser.add_argument("--max-contamination", type=float, default=10.0)
    parser.add_argument("--suffix", default=".fa")
    args = parser.parse_args(argv)

    kept, total = filter_bins(
        args.bins, args.report, args.outdir, args.out_report,
        args.max_contamination, args.suffix,
    )
    print(
        f"kept {kept} of {total} bins with contamination <= {args.max_contamination:g}%",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
