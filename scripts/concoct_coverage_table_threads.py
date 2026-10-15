"""
A script to generate the input coverage table for CONCOCT using a BEDFile.

Output is written to stdout. The BEDFile defines the regions used as
subcontigs for concoct. This makes it possible to get the coverage for
subcontigs without specifically mapping reads against the subcontigs.

@author: inodb, alneberg

Modified with Claude AI to add multithreading support.
"""

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO

# import glob
from signal import SIG_DFL, SIGPIPE, signal

import pandas as pd

# import tempfile


class BamProcessingError(Exception):
    """Raised when no BAM files could be successfully processed."""


class ContigMismatchError(Exception):
    """Raised when contig order differs between BAM files."""


class SampleNameCountMismatchError(Exception):
    """Raised when the number of sample names doesn't match the number of BAM files."""


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
                    f"ERROR! Something is wrong with the line:\n'{line}'\n"
                    "Perhaps 'concoct_part_' is misplaced or missing? Exiting!\n"
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

        fh = StringIO(out.decode("utf-8"))
        df = pd.read_csv(fh, header=None, sep="\t")

        # Calculate average coverage depth
        avg_coverage = df[4].divide(df[2] - df[1])

        return {
            "samplename": samplename,
            "contig_ids": df[3].values,
            "coverage": avg_coverage.values,
        }

    except (OSError, subprocess.SubprocessError, ValueError, pd.errors.ParserError) as e:
        sys.stderr.write(f"Exception processing {bamfile}: {e!s}\n")
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
            bamfile, _samplename = future_to_bam[future]
            try:
                result = future.result()
                if result is not None:
                    results.append(result)
                    sys.stderr.write(f"Completed processing: {bamfile}\n")
                else:
                    sys.stderr.write(f"Failed to process: {bamfile}\n")
            except (OSError, subprocess.SubprocessError, ValueError) as e:
                sys.stderr.write(f"Exception for {bamfile}: {e!s}\n")

    if not results:
        raise BamProcessingError("No BAM files were successfully processed")

    # Build the dataframe from results
    # Use the first result to establish contig order
    contig_ids = results[0]["contig_ids"]

    # Create dataframe with contigs as index
    df_dict = {}
    for result in results:
        # Ensure contig order matches
        if not all(result["contig_ids"] == contig_ids):
            raise ContigMismatchError("Contig order mismatch between BAM files")

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
        with open(args.samplenames) as sn_fh:
            samplenames = [s[:-1] for s in sn_fh]
        if len(samplenames) != len(args.bamfiles):
            raise SampleNameCountMismatchError(
                f"Nr of names ({len(samplenames)}) in samplenames should be equal to nr of given "
                f"bamfiles ({len(args.bamfiles)})"
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