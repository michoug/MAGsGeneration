"""Tests for the binner ranking/selection logic in scripts/rank_binners.py."""

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))

from rank_binners import (
    count_tiers,
    parse_quality_report_path,
    rank_score,
    select,
)

TIERS = {
    "mq": {"completeness": 50, "contamination": 10},
    "nc": {"completeness": 90, "contamination": 5},
    "nc_strict": {"completeness": 95, "contamination": 5},
}


def counts_row(sample, binner, mq, nc=0, nc_strict=0, n_bins=None):
    return {
        "sample": sample,
        "binner": binner,
        "n_bins": n_bins if n_bins is not None else mq,
        "mq": mq,
        "nc": nc,
        "nc_strict": nc_strict,
    }


class TestParseQualityReportPath:
    def test_extracts_sample_and_binner(self):
        path = "/x/Bins/intermediate/336R_mg_25/BinnerEval/comebin/checkm2/quality_report.tsv"
        sample, binner = parse_quality_report_path(path)
        assert sample == "336R_mg_25"
        assert binner == "comebin"

    def test_raises_on_unrecognized_path(self):
        with pytest.raises(ValueError):
            parse_quality_report_path("/not/a/matching/path.tsv")


class TestCountTiers:
    def test_counts_bins_per_tier(self):
        df = pd.DataFrame(
            {
                "Name": ["a", "b", "c"],
                "Completeness": [96, 92, 60],
                "Contamination": [1, 4, 9],
            }
        )
        counts = count_tiers(df, tiers=TIERS)
        assert counts["n_bins"] == 3
        assert counts["mq"] == 3
        assert counts["nc"] == 2
        assert counts["nc_strict"] == 1

    def test_empty_report_counts_zero(self):
        df = pd.DataFrame(columns=["Name", "Completeness", "Contamination"])
        counts = count_tiers(df, tiers=TIERS)
        assert counts == {"n_bins": 0, "mq": 0, "nc": 0, "nc_strict": 0}

    def test_completeness_thresholds_are_inclusive(self):
        # A bin landing exactly on a completeness cutoff must count toward
        # that tier (matches report/MAGs_summary.qmd's >= convention).
        df = pd.DataFrame(
            {
                "Name": ["at_mq", "at_nc", "at_nc_strict"],
                "Completeness": [50, 90, 95],
                "Contamination": [1, 1, 1],
            }
        )
        counts = count_tiers(df, tiers=TIERS)
        assert counts["mq"] == 3
        assert counts["nc"] == 2
        assert counts["nc_strict"] == 1


class TestRankScore:
    def test_lower_score_is_better(self):
        counts_df = pd.DataFrame(
            [
                counts_row("s1", "good", mq=10, nc=8, nc_strict=5),
                counts_row("s1", "bad", mq=3, nc=1, nc_strict=0),
                counts_row("s2", "good", mq=12, nc=9, nc_strict=6),
                counts_row("s2", "bad", mq=4, nc=2, nc_strict=1),
            ]
        )
        scores = rank_score(counts_df, tier_names=("mq", "nc", "nc_strict"))
        scores = scores.set_index("binner")
        assert scores.loc["good", "score"] < scores.loc["bad", "score"]
        assert scores.loc["good", "score"] == 1.0

    def test_ties_get_min_rank(self):
        counts_df = pd.DataFrame(
            [
                counts_row("s1", "a", mq=5, nc=5, nc_strict=5),
                counts_row("s1", "b", mq=5, nc=5, nc_strict=5),
                counts_row("s1", "c", mq=1, nc=1, nc_strict=1),
            ]
        )
        scores = rank_score(counts_df).set_index("binner")
        assert scores.loc["a", "score"] == scores.loc["b", "score"] == 1.0
        assert scores.loc["c", "score"] == 3.0


class TestSelect:
    CATALOGUE = {
        "classical_a": {"family": "classical"},
        "classical_b": {"family": "classical"},
        "ml_a": {"family": "ml_clustering"},
        "deep_a": {"family": "deep_representation"},
        "deep_b": {"family": "deep_representation"},
    }

    def test_one_per_family_prefers_best_of_each(self):
        counts_df = pd.DataFrame(
            [
                counts_row("s1", "classical_a", mq=20, nc=10),
                counts_row("s1", "classical_b", mq=15, nc=8),
                counts_row("s1", "ml_a", mq=18, nc=9),
                counts_row("s1", "deep_a", mq=25, nc=15),
                counts_row("s1", "deep_b", mq=22, nc=12),
            ]
        )
        result = select(counts_df, self.CATALOGUE, n=3, one_per_family=True)
        selected = result["selected"]
        assert len(selected) == 3
        assert "classical_b" not in selected  # classical_a is better in its family
        assert "deep_b" not in selected  # deep_a is better in its family
        assert set(selected) == {"classical_a", "ml_a", "deep_a"}

    def test_zero_mq_binner_is_excluded(self):
        counts_df = pd.DataFrame(
            [
                counts_row("s1", "classical_a", mq=20, nc=10),
                counts_row("s1", "ml_a", mq=0, nc=0),  # never produced an MQ bin
                counts_row("s1", "deep_a", mq=25, nc=15),
            ]
        )
        result = select(counts_df, self.CATALOGUE, n=3, one_per_family=True)
        assert "ml_a" not in result["selected"]
        report = result["report"].set_index("binner")
        assert report.loc["ml_a", "selected"] == False  # noqa: E712

    def test_backfill_when_fewer_families_than_n(self):
        # Only two families survive (one dropped for zero MQ); n=3 should
        # still return 3 binners by backfilling from the remaining pool.
        counts_df = pd.DataFrame(
            [
                counts_row("s1", "classical_a", mq=20, nc=10),
                counts_row("s1", "classical_b", mq=15, nc=8),
                counts_row("s1", "ml_a", mq=0, nc=0),
                counts_row("s1", "deep_a", mq=25, nc=15),
            ]
        )
        result = select(counts_df, self.CATALOGUE, n=3, one_per_family=True)
        assert len(result["selected"]) == 3
        assert set(result["selected"]) == {"classical_a", "classical_b", "deep_a"}

    def test_fewer_binners_than_n_returns_all_survivors(self):
        counts_df = pd.DataFrame(
            [
                counts_row("s1", "classical_a", mq=20, nc=10),
                counts_row("s1", "deep_a", mq=0, nc=0),
            ]
        )
        result = select(counts_df, self.CATALOGUE, n=3, one_per_family=True)
        assert result["selected"] == ["classical_a"]

    def test_pure_top_n_without_family_constraint(self):
        counts_df = pd.DataFrame(
            [
                counts_row("s1", "classical_a", mq=20, nc=10),
                counts_row("s1", "deep_a", mq=25, nc=15),
                counts_row("s1", "deep_b", mq=22, nc=12),
                counts_row("s1", "ml_a", mq=5, nc=1),
            ]
        )
        result = select(counts_df, self.CATALOGUE, n=3, one_per_family=False)
        assert set(result["selected"]) == {"deep_a", "deep_b", "classical_a"}
