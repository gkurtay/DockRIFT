from __future__ import annotations

import pandas as pd


def add_seed_ranks(x: pd.DataFrame) -> pd.DataFrame:
    y = x.copy()
    y["rank"] = y.groupby(["receptor", "seed"])["score"].rank(
        method="average", ascending=True
    )
    return y


def consensus_table(x: pd.DataFrame) -> pd.DataFrame:
    """Median-score consensus and per-ligand seed-rank envelope."""
    y = add_seed_ranks(x)
    s = (
        y.groupby(["ligand_id", "receptor"], as_index=False)
        .agg(
            score_median=("score", "median"),
            score_sd=("score", "std"),
            rank_min=("rank", "min"),
            rank_max=("rank", "max"),
        )
    )
    s["consensus_rank"] = s.groupby("receptor")["score_median"].rank(
        method="average", ascending=True
    )
    return s
