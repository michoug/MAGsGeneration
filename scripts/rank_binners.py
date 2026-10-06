#!/usr/bin/env python3
"""Rank binners on CheckM2 quality reports and select the best `n`.

Ranking follows Han, Wang & Zhu, Nat Commun 16:2865 (2025)
(https://doi.org/10.1038/s41467-025-57957-6): for each sample, binners are
ranked (1 = best) on the count of bins meeting each quality tier; a binner's
score is the mean of its per-tier mean rank across samples (lower is
better). The paper's HQ tier (rRNA + tRNA presence) is approximated here by
a stricter completeness/contamination tier (`nc_strict`) to avoid adding
barrnap/tRNAscan-SE to the pipeline.
"""

import argparse
import re
import sys

import pandas as pd

DEFAULT_TIERS = {
    "mq": {"completeness": 50, "contamination": 10},
    "nc": {"completeness": 90, "contamination": 5},
    "nc_strict": {"completeness": 95, "contamination": 5},
}

PATH_RE = re.compile(r"(?P<sample>[^/]+)/BinnerEval/(?P<binner>[^/]+)/checkm2/")


def parse_quality_report_path(path):
    match = PATH_RE.search(path.replace("\\", "/"))
    if not match:
        raise ValueError(f"cannot infer sample/binner from path: {path}")
    return match.group("sample"), match.group("binner")


def count_tiers(quality_df, tiers=None):
    tiers = tiers or DEFAULT_TIERS
    counts = {"n_bins": len(quality_df)}
    for name, thresh in tiers.items():
        if len(quality_df) == 0:
            counts[name] = 0
            continue
        mask = (quality_df["Completeness"] > thresh["completeness"]) & (
            quality_df["Contamination"] < thresh["contamination"]
        )
        counts[name] = int(mask.sum())
    return counts


def load_quality_report(path):
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["Name", "Completeness", "Contamination"])


def load_counts(paths, tiers=None):
    rows = []
    for path in paths:
        sample, binner = parse_quality_report_path(path)
        counts = count_tiers(load_quality_report(path), tiers=tiers)
        rows.append({"sample": sample, "binner": binner, **counts})
    return pd.DataFrame(rows)


def rank_score(counts_df, tier_names=("mq", "nc", "nc_strict")):
    """Lower score is better. One row per binner."""
    df = counts_df.copy()
    rank_cols = []
    for tier in tier_names:
        col = f"rank_{tier}"
        df[col] = df.groupby("sample")[tier].rank(method="min", ascending=False)
        rank_cols.append(col)

    per_binner = df.groupby("binner", as_index=False)[rank_cols].mean()
    per_binner["score"] = per_binner[rank_cols].mean(axis=1)
    totals = counts_df.groupby("binner", as_index=False)[
        ["n_bins", *tier_names]
    ].sum()
    out = per_binner.merge(totals, on="binner")
    return out.sort_values("score").reset_index(drop=True)


def select(counts_df, catalogue, n=3, one_per_family=True, tier_names=("mq", "nc", "nc_strict")):
    """Returns {"selected": [binner, ...], "report": DataFrame}."""
    all_binners = sorted(counts_df["binner"].unique())
    mq_totals = counts_df.groupby("binner")["mq"].sum()
    survivors = sorted(mq_totals[mq_totals > 0].index.tolist())
    dropped = sorted(set(all_binners) - set(survivors))

    def family_of(binner):
        return catalogue.get(binner, {}).get("family", "unknown")

    if survivors:
        scored = rank_score(counts_df[counts_df["binner"].isin(survivors)], tier_names=tier_names)
    else:
        scored = pd.DataFrame(columns=["binner", "score", "n_bins", *tier_names])
    scored["family"] = scored["binner"].map(family_of)

    selected = []
    if one_per_family:
        best_per_family = scored.sort_values("score").drop_duplicates(
            "family", keep="first"
        )
        selected = best_per_family.sort_values("score")["binner"].tolist()[:n]
        if len(selected) < n:
            remaining = scored[~scored["binner"].isin(selected)].sort_values("score")
            selected += remaining["binner"].tolist()[: n - len(selected)]
    else:
        selected = scored.sort_values("score")["binner"].tolist()[:n]

    report_rows = scored.to_dict("records")
    for binner in dropped:
        report_rows.append(
            {
                "binner": binner,
                "family": family_of(binner),
                "score": float("nan"),
                "n_bins": int(mq_totals.get(binner, 0)) if binner in mq_totals else 0,
                **{t: 0 for t in tier_names},
            }
        )
    report = pd.DataFrame(report_rows)
    report["selected"] = report["binner"].isin(selected)
    report = report.sort_values(["selected", "score"], ascending=[False, True])

    return {"selected": selected, "report": report}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quality-reports", nargs="+", required=True)
    parser.add_argument("--catalogue-json", required=True, help="JSON mapping binner -> {family, label}")
    parser.add_argument("--n", type=int, default=3)
    parser.add_argument("--one-per-family", action="store_true", default=True)
    parser.add_argument("--no-one-per-family", dest="one_per_family", action="store_false")
    parser.add_argument("--out-ranking", required=True)
    parser.add_argument("--out-report", required=True)
    parser.add_argument("--out-selected", required=True)
    args = parser.parse_args(argv)

    import json

    with open(args.catalogue_json) as fh:
        catalogue = json.load(fh)

    counts_df = load_counts(args.quality_reports)
    counts_df.to_csv(args.out_ranking, sep="\t", index=False)

    result = select(counts_df, catalogue, n=args.n, one_per_family=args.one_per_family)
    result["report"].to_csv(args.out_report, sep="\t", index=False)

    with open(args.out_selected, "w") as fh:
        fh.writelines(binner + "\n" for binner in result["selected"])

    for binner in result["selected"]:
        row = result["report"][result["report"]["binner"] == binner].iloc[0]
        print(
            f"selected {binner}: family={row['family']} score={row['score']:.3f}",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
