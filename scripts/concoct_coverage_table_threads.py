#!/usr/bin/env python
"""
A script to generate the input coverage table for CONCOCT using a BEDFile.

Output is written to stdout. The BEDFile defines the regions used as
subcontigs for concoct. This makes it possible to get the coverage for
subcontigs without specifically mapping reads against the subcontigs.

@author: inodb, alneberg

Modified with Claude AI to add multithreading support.
"""

from __future__ import print_function
import sys
import os
import argparse
import subprocess
# import glob
from signal import signal, SIGPIPE, SIG_DFL
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
# import tempfile


def check_bed_file_for_errors(bedfile):
    with open(bedfile) as ifh:
        for line in ifh:
            line = line.strip()
            original_id, _, _, cutup_id = line.split("\t")
            try:
                assert "concoct_part_" not in original_id
                assert "concoct_part_" in cutup_id
            except AssertionError:
                sys.stderr.write(
                    (
                        "ERROR! Something is wrong with the line:\n'{}'\n"
                        "Perhaps 'concoct_part_' is misplaced or missing? Exiting!\n"
                    ).format(line)
                )
                sys.exit(-1)


def process_single_bam(bedfile, bamfile, samplename):
    """Process a single BAM file and return its coverage data."""
    try:
        p = subprocess.Popen(
            ["samtools", "bedcov", bedfile, bamfile],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        out, err = p.communicate()

        if p.returncode != 0:
            sys.stderr.write(
                "Error processing {}: {}\n".format(bamfile, err.decode("utf-8"))
            )
            return None

        # Parse output
        if sys.version_info[0] < 3:
            from StringIO import StringIO
        else:
            from io import StringIO

        fh = StringIO(out.decode("utf-8"))
        df = pd.read_csv(fh, header=None, sep="\t")

        # Calculate average coverage depth
        avg_coverage = df[4].divide(df[2] - df[1])

        return {
            "samplename": samplename,
            "contig_ids": df[3].values,
            "coverage": avg_coverage.values,
        }

    except Exception as e:
        sys.stderr.write("Exception processing {}: {}\n".format(bamfile, str(e)))
        return None


def generate_input_table(bedfile, bamfiles, samplenames=None, threads=1):
    """Reads input files into dictionaries then prints everything in the table
    format required for running CONCOCT."""

    # Prepare sample names
    if samplenames is None:
        col_names = [
            os.path.splitext(os.path.basename(bamfile))[0] for bamfile in bamfiles
        ]
    else:
        col_names = samplenames

    # Process BAM files in parallel
    results = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        # Submit all tasks
        future_to_bam = {
            executor.submit(process_single_bam, bedfile, bamfile, samplename): (
                bamfile,
                samplename,
            )
            for bamfile, samplename in zip(bamfiles, col_names)
        }

        # Collect results as they complete
        for future in as_completed(future_to_bam):
            bamfile, samplename = future_to_bam[future]
            try:
                result = future.result()
                if result is not None:
                    results.append(result)
                    sys.stderr.write("Completed processing: {}\n".format(bamfile))
                else:
                    sys.stderr.write("Failed to process: {}\n".format(bamfile))
            except Exception as e:
                sys.stderr.write("Exception for {}: {}\n".format(bamfile, str(e)))

    if not results:
        raise Exception("No BAM files were successfully processed")

    # Build the dataframe from results
    # Use the first result to establish contig order
    contig_ids = results[0]["contig_ids"]

    # Create dataframe with contigs as index
    df_dict = {}
    for result in results:
        # Ensure contig order matches
        if not all(result["contig_ids"] == contig_ids):
            raise Exception("Contig order mismatch between BAM files")

        col_name = "cov_mean_sample_{}".format(result["samplename"])
        df_dict[col_name] = result["coverage"]

    df = pd.DataFrame(df_dict, index=contig_ids)
    df.index.name = "contig"

    # Output to stdout
    df.to_csv(sys.stdout, sep="\t", float_format="%.3f")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "bedfile",
        help="Contigs BEDFile with four columns representing: 'Contig ID, Start Position, "
        "End Position and SubContig ID' respectively. The Subcontig ID must contain the pattern "
        "'concoct_part_[0-9]*' while the contigs which are not cutup cannot contain this pattern. "
        "This file can be generated by the cut_up_fasta.py script.",
    )
    parser.add_argument(
        "bamfiles", nargs="+", help="BAM files with mappings to the original contigs."
    )
    parser.add_argument(
        "--samplenames",
        default=None,
        help="File with sample names, one line each. Should be same nr "
        "of bamfiles. Default sample names used are the file names of the bamfiles, excluding the file extension.",
    )
    parser.add_argument(
        "--threads",
        "-t",
        type=int,
        default=1,
        help="Number of threads to use for parallel processing (default: 1)",
    )
    args = parser.parse_args()

    # Get sample names
    if args.samplenames is not None:
        samplenames = [s[:-1] for s in open(args.samplenames).readlines()]
        if len(samplenames) != len(args.bamfiles):
            raise Exception(
                "Nr of names ({0}) in samplenames should be equal to nr of given "
                "bamfiles ({1})".format(len(samplenames), len(args.bamfiles))
            )
    else:
        samplenames = None

    # ignore broken pipe error when piping output
    # http://newbebweb.blogspot.pt/2012/02/python-head-ioerror-errno-32-broken.html
    signal(SIGPIPE, SIG_DFL)
    check_bed_file_for_errors(args.bedfile)
    generate_input_table(
        args.bedfile, args.bamfiles, samplenames=samplenames, threads=args.threads
    )
