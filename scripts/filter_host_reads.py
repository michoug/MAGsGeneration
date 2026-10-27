"""Drop read pairs that genuinely align to a host genome, keep everything else.

Reads a name-grouped SAM stream (strobealign output, mates adjacent) on stdin
and writes the SAM records of the non-host pairs to stdout.

`samtools view -f 12` (keep pairs with both mates unmapped) is not usable on
strobealign output: strobealign reports local hits of ~20 bp with the rest of
the read soft-clipped, often with MAPQ 60 (strobealign issue #477). Against
hg38 this flagged ~50% of river metagenome pairs as human while only ~0.1%
of reads actually aligned to it. A mate is therefore counted as host only if
its primary alignment covers at least `--min-aligned-fraction` of the read
with at most `--max-mismatch-rate` edit distance (NM) per aligned base; a pair
is dropped if either mate is host.
"""

import argparse
import re
import sys

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
UNMAPPED, SECONDARY, SUPPLEMENTARY = 0x4, 0x100, 0x800


def is_host(fields, min_aligned_fraction, max_mismatch_rate):
    """True if this primary SAM record is a real alignment to the host."""
    if int(fields[1]) & UNMAPPED:
        return False
    read_len = len(fields[9])
    if read_len == 0 or fields[9] == "*":
        return False
    aligned = sum(int(n) for n, op in CIGAR_RE.findall(fields[5]) if op in "M=X")
    if aligned < min_aligned_fraction * read_len:
        return False
    nm = next((int(t[5:]) for t in fields[11:] if t.startswith("NM:i:")), 0)
    return nm <= max_mismatch_rate * aligned


def filter_sam(lines, out, min_aligned_fraction=0.8, max_mismatch_rate=0.05):
    """Write header and non-host records to *out*. Returns (kept, dropped) pairs."""
    kept = dropped = 0
    group, group_name, host = [], None, False

    def flush():
        nonlocal kept, dropped
        if not group:
            return
        if host:
            dropped += 1
        else:
            kept += 1
            out.writelines(group)

    for line in lines:
        if line.startswith("@"):
            out.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if int(fields[1]) & (SECONDARY | SUPPLEMENTARY):
            continue
        if fields[0] != group_name:
            flush()
            group, group_name, host = [], fields[0], False
        group.append(line)
        host = host or is_host(fields, min_aligned_fraction, max_mismatch_rate)
    flush()
    return kept, dropped


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-aligned-fraction", type=float, default=0.8)
    parser.add_argument("--max-mismatch-rate", type=float, default=0.05)
    args = parser.parse_args(argv)

    kept, dropped = filter_sam(
        sys.stdin, sys.stdout, args.min_aligned_fraction, args.max_mismatch_rate
    )
    total = kept + dropped
    print(
        f"kept {kept} pairs, removed {dropped} host pairs"
        f" ({100 * dropped / total if total else 0:.3f}%)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
