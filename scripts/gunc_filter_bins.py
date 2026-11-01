"""Keep refined bins that pass GUNC and a contamination ceiling.

Reads GUNC's maxCSS_level table and a CheckM2-format quality report
(Name / Completeness / Contamination) for the bins in --bins, copies the
bins with pass.GUNC not False (True, or nan when too few genes map for GUNC
to decide) and Contamination <= --max-contamination to --outdir, and writes
one row per bin with the decision.
"""

import argparse
import csv
import os
import shutil
import sys

COLUMNS = [
    "bin",
    "completeness",
    "contamination",
    "gunc_css",
    "gunc_contamination_portion",
    "pass_gunc",
    "kept",
    "reason",
]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def decide(quality, gunc, max_contamination):
    """(kept, reason) for one bin from its quality row and GUNC row (or None)."""
    reasons = []
    if gunc is not None and gunc.get("pass.GUNC") == "False":
        reasons.append("gunc_chimeric")
    if float(quality["Contamination"]) > max_contamination:
        reasons.append(f"contamination>{max_contamination:g}")
    return not reasons, ";".join(reasons) or "pass"


def filter_bins(bins_dir, quality_path, gunc_path, outdir, table_path, suffix=".fa",
                max_contamination=10.0):
    """Copy kept bins to *outdir*, write the decision table; return (kept, total)."""
    quality = {r["Name"]: r for r in read_tsv(quality_path)}
    gunc = {r["genome"]: r for r in read_tsv(gunc_path)}
    os.makedirs(outdir, exist_ok=True)
    kept_n = total = 0
    with open(table_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for fn in sorted(os.listdir(bins_dir)):
            if not fn.endswith(suffix):
                continue
            name = fn[: -len(suffix)]
            q = quality[name]
            g = gunc.get(name)
            kept, reason = decide(q, g, max_contamination)
            total += 1
            if kept:
                kept_n += 1
                shutil.copy(os.path.join(bins_dir, fn), os.path.join(outdir, fn))
            w.writerow([
                name,
                q["Completeness"],
                q["Contamination"],
                g.get("clade_separation_score", "") if g else "",
                g.get("contamination_portion", "") if g else "",
                g.get("pass.GUNC", "") if g else "",
                "yes" if kept else "no",
                reason,
            ])
    return kept_n, total


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bins", required=True)
    parser.add_argument("--quality", required=True, help="CheckM2-format report")
    parser.add_argument("--gunc", required=True, help="GUNC maxCSS_level.tsv")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--table", required=True)
    parser.add_argument("--suffix", default=".fa")
    parser.add_argument("--max-contamination", type=float, default=10.0)
    args = parser.parse_args(argv)

    kept, total = filter_bins(
        args.bins, args.quality, args.gunc, args.outdir, args.table,
        args.suffix, args.max_contamination,
    )
    print(f"kept {kept} of {total} bins", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
