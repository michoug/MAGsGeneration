#!/usr/bin/env python3
"""Split a DAS-Tool-style contig2bin TSV into per-bin FASTA files.

Used by the binner_selection step to re-materialize each binner's raw bins
(from the contig2bin tables already produced by prepare_dasTool) so they can
be scored independently with CheckM2.
"""

import argparse
import os
import sys


def parse_fasta(path):
    sequences = {}
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            sequences[name] = "".join(chunks)
    return sequences


def read_contig2bin(path):
    bins = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            contig, bin_id = parts[0], parts[1]
            bins.setdefault(bin_id, []).append(contig)
    return bins


def write_fasta(path, records, width=80):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            fh.writelines(seq[i : i + width] + "\n" for i in range(0, len(seq), width))


def explode(contig2bin_path, assembly_path, outdir, min_size=0):
    """Write one FASTA per bin; returns (n_bins_written, n_contigs_missing)."""
    os.makedirs(outdir, exist_ok=True)

    bins = read_contig2bin(contig2bin_path)
    if not bins:
        return 0, 0

    sequences = parse_fasta(assembly_path)

    n_written = 0
    n_missing = 0
    for bin_id, contigs in bins.items():
        records = []
        for contig in contigs:
            seq = sequences.get(contig)
            if seq is None:
                n_missing += 1
                continue
            records.append((contig, seq))
        if sum(len(seq) for _, seq in records) < min_size:
            continue
        write_fasta(os.path.join(outdir, f"{bin_id}.fa"), records)
        n_written += 1

    return n_written, n_missing


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contig2bin", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-size", type=int, default=0)
    args = parser.parse_args(argv)

    n_written, n_missing = explode(
        args.contig2bin, args.assembly, args.outdir, args.min_size
    )
    if n_missing:
        print(
            f"error: {n_missing} contig IDs in {args.contig2bin} "
            f"were not found in {args.assembly} - contig2bin table and "
            "assembly are out of sync, refusing to score a truncated bin",
            file=sys.stderr,
        )
        return 1
    print(f"wrote {n_written} bin(s) to {args.outdir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
