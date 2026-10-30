"""Tests for scripts/filter_host_reads.py."""

import io
import os
import sys

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"),
)

from filter_host_reads import filter_sam

SEQ = "A" * 150
HEADER = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000000\n"


def rec(name, flag, cigar="*", nm=None, pos=0, rname="*"):
    tags = f"\tNM:i:{nm}" if nm is not None else ""
    return f"{name}\t{flag}\t{rname}\t{pos}\t60\t{cigar}\t*\t0\t0\t{SEQ}\t{'I' * 150}{tags}\n"


def run(records, **kw):
    out = io.StringIO()
    kept, dropped = filter_sam(io.StringIO(HEADER + "".join(records)), out, **kw)
    body = [l for l in out.getvalue().splitlines() if not l.startswith("@")]
    return kept, dropped, body, out.getvalue()


class TestFilterSam:
    def test_keeps_header(self):
        _, _, _, text = run([rec("r1", 77), rec("r1", 141)])
        assert text.startswith(HEADER)

    def test_both_unmapped_pair_is_kept(self):
        kept, dropped, body, _ = run([rec("r1", 77), rec("r1", 141)])
        assert (kept, dropped) == (1, 0)
        assert len(body) == 2

    def test_short_softclipped_hit_is_not_host(self):
        # the strobealign artefact: 21 bp aligned, rest soft-clipped, MAPQ 60
        kept, dropped, body, _ = run(
            [rec("r1", 73, "53S21M76S", 0, 100, "chr1"), rec("r1", 133)]
        )
        assert (kept, dropped) == (1, 0)
        assert len(body) == 2

    def test_full_length_alignment_removes_the_pair(self):
        kept, dropped, body, _ = run(
            [rec("r1", 73, "150M", 2, 100, "chr1"), rec("r1", 133)]
        )
        assert (kept, dropped) == (0, 1)
        assert body == []

    def test_full_length_but_divergent_is_not_host(self):
        kept, dropped, _, _ = run([rec("r1", 73, "150M", 20, 100, "chr1"), rec("r1", 133)])
        assert (kept, dropped) == (1, 0)

    def test_secondary_records_are_dropped_and_ignored(self):
        kept, dropped, body, _ = run(
            [rec("r1", 77), rec("r1", 141), rec("r1", 329, "150M", 0, 5, "chr1")]
        )
        assert (kept, dropped) == (1, 0)
        assert len(body) == 2

    def test_pairs_are_judged_independently(self):
        kept, dropped, body, _ = run(
            [
                rec("a", 77), rec("a", 141),
                rec("b", 67, "150M", 0, 10, "chr1"), rec("b", 131, "150M", 0, 300, "chr1"),
                rec("c", 77), rec("c", 141),
            ]
        )
        assert (kept, dropped) == (2, 1)
        assert [l.split("\t")[0] for l in body] == ["a", "a", "c", "c"]
