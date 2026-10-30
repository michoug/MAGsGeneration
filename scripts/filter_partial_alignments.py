"""Drop short partial alignments from a SAM stream before coverage is computed.

strobealign reports local hits where most of the read is soft-clipped
(strobealign issue #477). When reads are mapped to an assembly these add
spurious coverage: on the Nant BAMs, 21-36% of mapped reads aligned over less
than 80% of their length in the middle of a contig (median 40-57 bp).

A mapped record is kept if its alignment covers at least
`--min-aligned-fraction` of the read, or if it starts or ends within
`--edge` bp of a contig end (reads overhanging a contig end are legitimately
partial). Unmapped records and the header pass through unchanged; secondary
and supplementary records follow the same rule.
"""

import argparse
import sys

from filter_host_reads import CIGAR_RE, UNMAPPED


def keep_record(fields, contig_len, min_aligned_fraction, edge):
    if int(fields[1]) & UNMAPPED:
        return True
    ops = CIGAR_RE.findall(fields[5])
    aligned = sum(int(n) for n, op in ops if op in "M=X")
    read_len = sum(int(n) for n, op in ops if op in "MIS=X")
    if read_len == 0 or aligned >= min_aligned_fraction * read_len:
        return True
    start = int(fields[3])
    end = start + sum(int(n) for n, op in ops if op in "MDN=X") - 1
    length = contig_len.get(fields[2])
    return start <= 1 + edge or (length is not None and end >= length - edge)


def filter_sam(lines, out, min_aligned_fraction=0.8, edge=5):
    """Write header and kept records to *out*. Returns (kept, dropped) records."""
    contig_len = {}
    kept = dropped = 0
    for line in lines:
        if line.startswith("@"):
            if line.startswith("@SQ"):
                tags = dict(t.split(":", 1) for t in line.rstrip("\n").split("\t")[1:])
                contig_len[tags["SN"]] = int(tags["LN"])
            out.write(line)
            continue
        fields = line.split("\t", 11)
        if keep_record(fields, contig_len, min_aligned_fraction, edge):
            kept += 1
            out.write(line)
        else:
            dropped += 1
    return kept, dropped


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-aligned-fraction", type=float, default=0.8)
    parser.add_argument("--edge", type=int, default=5)
    args = parser.parse_args(argv)

    kept, dropped = filter_sam(sys.stdin, sys.stdout, args.min_aligned_fraction, args.edge)
    total = kept + dropped
    print(
        f"kept {kept} records, removed {dropped} partial alignments"
        f" ({100 * dropped / total if total else 0:.2f}%)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
