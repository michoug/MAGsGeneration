"""Tests for scripts/filter_partial_alignments.py."""

import io
import os
import sys

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"),
)

from filter_partial_alignments import filter_sam

HEADER = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:2000\n"
SEQ = "A" * 150


def rec(flag, cigar, pos, rname="c1"):
    return f"r\t{flag}\t{rname}\t{pos}\t60\t{cigar}\t*\t0\t0\t{SEQ}\t*\n"


def run(records):
    out = io.StringIO()
    kept, dropped = filter_sam(io.StringIO(HEADER + "".join(records)), out)
    return kept, dropped, out.getvalue()


class TestFilterSam:
    def test_header_passes_through(self):
        _, _, text = run([])
        assert text == HEADER

    def test_full_length_alignment_is_kept(self):
        assert run([rec(0, "150M", 500)])[:2] == (1, 0)

    def test_interior_partial_hit_is_dropped(self):
        kept, dropped, text = run([rec(0, "53S21M76S", 500)])
        assert (kept, dropped) == (0, 1)
        assert text == HEADER

    def test_overhang_at_contig_start_is_kept(self):
        assert run([rec(0, "60S90M", 1)])[:2] == (1, 0)

    def test_overhang_at_contig_end_is_kept(self):
        # 1911 + 90 - 1 = 2000, the last base of c1
        assert run([rec(0, "90M60S", 1911)])[:2] == (1, 0)

    def test_unmapped_record_is_kept(self):
        assert run([rec(4, "*", 0, "*")])[:2] == (1, 0)

    def test_small_indels_do_not_count_as_partial(self):
        assert run([rec(0, "70M2I78M", 500)])[:2] == (1, 0)
