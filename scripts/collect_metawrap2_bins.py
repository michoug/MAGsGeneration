"""Collect metawrap2's refined bins under stable per-sample names, dropping GUNC chimeras.

metawrap2 names consolidated bins bin_001.fasta, bin_002.fasta, ... but keeps
the input names when only one bin set was given, and the same bin_001 exists
in every sample. Each bin listed in metawrap2's .stats.tsv is copied to
<outdir>/<sample>_bin_NNN.fa (numbered in .stats.tsv order).

metawrap2 only annotates GUNC's verdict (gunc_pass yes/no, empty when GUNC
could not assess a bin); bins with gunc_pass "no" are not copied. Every bin is
listed in the --table output with its original name, new name, quality and
GUNC verdict.
"""

import argparse
import csv
import os
import shutil
import sys

BIN_EXTENSION = ".fasta"
OUT_EXTENSION = ".fa"
TABLE_COLUMNS = [
    "new_name",
    "metawrap2_bin",
    "completeness",
    "contamination",
    "gunc_css",
    "gunc_pass",
    "kept",
]


def read_stats(path):
    """Rows of metawrap2's .stats.tsv as dicts; [] if it is missing or empty."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def collect(bins_dir, stats_path, outdir, table_path, sample, drop_chimeric=True):
    """Copy the kept bins and write the name table. Returns (kept, dropped)."""
    os.makedirs(outdir, exist_ok=True)
    rows = read_stats(stats_path)
    kept = dropped = 0
    with open(table_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(TABLE_COLUMNS)
        for row in rows:
            name = row["bin"]
            src = os.path.join(bins_dir, name + BIN_EXTENSION)
            if not os.path.isfile(src):
                raise FileNotFoundError(f"bin listed in {stats_path} is missing: {src}")
            gunc_pass = (row.get("gunc_pass") or "").strip().lower()
            chimeric = drop_chimeric and gunc_pass == "no"
            new_name = ""
            if chimeric:
                dropped += 1
            else:
                kept += 1
                new_name = f"{sample}_bin_{kept:03d}"
                shutil.copy(src, os.path.join(outdir, new_name + OUT_EXTENSION))
            writer.writerow(
                [
                    new_name,
                    name,
                    row.get("completeness", ""),
                    row.get("contamination", ""),
                    row.get("gunc_css", ""),
                    gunc_pass,
                    "no" if chimeric else "yes",
                ]
            )
    return kept, dropped


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bins", required=True, help="metawrap_<c>_<x>_bins folder")
    parser.add_argument("--stats", required=True, help="metawrap_<c>_<x>_bins.stats.tsv")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--table", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--keep-chimeric",
        action="store_true",
        help="copy bins that failed GUNC as well (they are still flagged in --table)",
    )
    args = parser.parse_args(argv)

    kept, dropped = collect(
        args.bins,
        args.stats,
        args.outdir,
        args.table,
        args.sample,
        drop_chimeric=not args.keep_chimeric,
    )
    print(
        f"{args.sample}: kept {kept} bins, dropped {dropped} that failed GUNC",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
