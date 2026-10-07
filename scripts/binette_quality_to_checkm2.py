#!/usr/bin/env python3
"""Reformat Binette's final_bins_quality_reports.tsv into CheckM2's native schema.

Binette runs CheckM2 internally to score and select its final bins, so
re-running `checkm2 predict` on Binette's output would redo the same work for
no new information. This script instead maps Binette's own report onto
CheckM2's native column layout (Name, Completeness, Contamination, ...) so
downstream rules (rosella refine --checkm-results, checkm2_refine, galah)
that expect that schema keep working unchanged. Only Name/Completeness/
Contamination are populated from real data; columns Binette doesn't compute
are filled with "NA".
"""

import argparse
import csv
import sys

CHECKM2_HEADER = [
    "Name",
    "Completeness",
    "Contamination",
    "Completeness_Model_Used",
    "Translation_Table_Used",
    "Coding_Density",
    "Contig_N50",
    "Average_Gene_Length",
    "Genome_Size",
    "GC_Content",
    "Total_Coding_Sequences",
    "Total_Contigs",
    "Max_Contig_Length",
    "Additional_Notes",
]

# CheckM2 column -> Binette column (omitted entries are filled with "NA")
FROM_BINETTE = {
    "Name": "name",
    "Completeness": "completeness",
    "Contamination": "contamination",
    "Completeness_Model_Used": "checkm2_model",
    "Coding_Density": "coding_density",
    "Contig_N50": "N50",
    "Genome_Size": "size",
    "Total_Contigs": "contig_count",
}


def convert(binette_rows):
    rows = []
    for row in binette_rows:
        out = {}
        for col in CHECKM2_HEADER:
            src = FROM_BINETTE.get(col)
            out[col] = row[src] if src is not None else "NA"
        rows.append(out)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binette-tsv", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    with open(args.binette_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = convert(reader)

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CHECKM2_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {len(rows)} row(s) to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
