from __future__ import annotations

import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

from .rrl import conservative_pair_rrl


def screening_metrics(labels, scores, *, topk: int = 20, lower_is_better: bool = True) -> dict:
    """Compute simple experimental-label screening metrics."""
    y = np.asarray(labels, int)
    s = np.asarray(scores, float)
    if y.ndim != 1 or s.ndim != 1 or len(y) != len(s):
        raise ValueError("labels and scores must be equal-length 1D arrays")
    if not set(np.unique(y)).issubset({0, 1}) or len(np.unique(y)) < 2:
        raise ValueError("labels must contain both binary classes 0 and 1")
    if not 1 <= topk <= len(y):
        raise ValueError("topk must lie between 1 and the number of ligands")

    pred = -s if lower_is_better else s
    order = np.argsort(s if lower_is_better else -s, kind="mergesort")
    top = order[:topk]
    active_total = int(y.sum())
    top_active = int(y[top].sum())
    expected = topk * active_total / len(y)
    ef = float(top_active / expected) if expected > 0 else np.nan
    return {
        "roc_auc": float(roc_auc_score(y, pred)),
        "pr_auc": float(average_precision_score(y, pred)),
        "topk": int(topk),
        "topk_active": top_active,
        "topk_active_fraction": float(top_active / topk),
        "ef_topk": ef,
    }


def pair_validity_instability(labels, scores_A, scores_B, *, topk: int = 20) -> dict:
    """Compare experimental-label screening performance across two receptor conformations."""
    a = screening_metrics(labels, scores_A, topk=topk)
    b = screening_metrics(labels, scores_B, topk=topk)
    y = np.asarray(labels, int)
    A = np.asarray(scores_A, float)
    B = np.asarray(scores_B, float)
    oa = np.argsort(A, kind="mergesort")[:topk]
    ob = np.argsort(B, kind="mergesort")[:topk]
    aa = set(np.asarray(oa)[y[oa] == 1].tolist())
    bb = set(np.asarray(ob)[y[ob] == 1].tolist())
    union = aa | bb
    jaccard = float(len(aa & bb) / len(union)) if union else np.nan
    return {
        "roc_auc_A": a["roc_auc"],
        "roc_auc_B": b["roc_auc"],
        "abs_delta_roc_auc": abs(a["roc_auc"] - b["roc_auc"]),
        "pr_auc_A": a["pr_auc"],
        "pr_auc_B": b["pr_auc"],
        "abs_delta_pr_auc": abs(a["pr_auc"] - b["pr_auc"]),
        "ef_topk_A": a["ef_topk"],
        "ef_topk_B": b["ef_topk"],
        "abs_delta_ef_topk": abs(a["ef_topk"] - b["ef_topk"]),
        "active_topk_jaccard": jaccard,
    }


def _reversal_rate_above(base, pert, threshold):
    base = np.asarray(base, float)
    pert = np.asarray(pert, float)
    ii, jj = np.triu_indices(len(base), 1)
    d0 = base[ii] - base[jj]
    d1 = pert[ii] - pert[jj]
    valid = (d0 != 0) & (d1 != 0)
    gap = np.abs(d0[valid])
    rev = np.sign(d0[valid]) != np.sign(d1[valid])
    sel = gap >= threshold
    n = int(sel.sum())
    return n, (float(rev[sel].mean()) if n else np.nan)


def heldout_rrl_calibration(scores_A, scores_B, *, repeats: int = 100,
                            train_fraction: float = 0.8, alpha: float = 0.05,
                            target_reversal: float | None = None,
                            rng_seed: int = 20260817) -> list[dict]:
    """Repeated ligand-identity train/test stress test for an observed-support RRL.

    The function deliberately treats calibration as a *diagnostic*. A finite training
    threshold can be non-evaluable in a held-out split if the test subset contains no
    score gaps at or above the learned threshold. Results should therefore not be
    interpreted as a universal predictive certificate.
    """
    A = np.asarray(scores_A, float)
    B = np.asarray(scores_B, float)
    if A.ndim != 1 or B.ndim != 1 or len(A) != len(B):
        raise ValueError("score arrays must be equal-length 1D arrays")
    if not (0 < train_fraction < 1):
        raise ValueError("train_fraction must lie in (0,1)")
    if repeats < 1:
        raise ValueError("repeats must be positive")
    if target_reversal is None:
        target_reversal = alpha
    n = len(A)
    n_train = int(round(train_fraction * n))
    if n_train < 3 or n - n_train < 3:
        raise ValueError("train/test split leaves too few ligands")
    rng = np.random.default_rng(rng_seed)
    rows = []
    for rep in range(repeats):
        p = rng.permutation(n)
        tr, te = p[:n_train], p[n_train:]
        r = conservative_pair_rrl(A[tr], B[tr], alpha=alpha)
        if r.conservative_status == "FINITE":
            thr = float(r.conservative)
            nA, rateA = _reversal_rate_above(A[te], B[te], thr)
            nB, rateB = _reversal_rate_above(B[te], A[te], thr)
            evaluable = nA > 0 and nB > 0
            passed = bool(rateA <= target_reversal and rateB <= target_reversal) if evaluable else None
            mx = max(rateA, rateB) if evaluable else np.nan
        else:
            thr = np.nan
            nA = nB = 0
            rateA = rateB = mx = np.nan
            evaluable = False
            passed = None
        rows.append({
            "repeat": rep + 1,
            "train_status": r.conservative_status,
            "train_rrl_conservative": thr,
            "test_A_pairs_above": nA,
            "test_B_pairs_above": nB,
            "test_A_reversal_rate": rateA,
            "test_B_reversal_rate": rateB,
            "max_test_reversal_rate": mx,
            "evaluable": evaluable,
            "passes_target_both_directions": passed,
        })
    return rows


def subsample_rrl_sensitivity(scores_A, scores_B, *, sample_sizes=(100, 200, 300),
                              repeats: int = 100, alpha: float = 0.05,
                              rng_seed: int = 20260818) -> list[dict]:
    """Quantify ligand-universe sensitivity by repeated without-replacement subsampling."""
    A = np.asarray(scores_A, float)
    B = np.asarray(scores_B, float)
    if A.ndim != 1 or B.ndim != 1 or len(A) != len(B):
        raise ValueError("score arrays must be equal-length 1D arrays")
    rng = np.random.default_rng(rng_seed)
    rows = []
    for size in sample_sizes:
        if not 3 <= int(size) <= len(A):
            raise ValueError("sample size outside valid range")
        for rep in range(repeats):
            ix = rng.choice(len(A), int(size), replace=False)
            r = conservative_pair_rrl(A[ix], B[ix], alpha=alpha)
            value = r.conservative if r.conservative_status == "FINITE" else r.conservative_lower_bound
            rows.append({
                "sample_size": int(size),
                "repeat": rep + 1,
                "status": r.conservative_status,
                "rrl_value_or_lower_bound": float(value),
            })
    return rows
