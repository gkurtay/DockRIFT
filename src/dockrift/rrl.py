from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from sklearn.isotonic import IsotonicRegression


@dataclass(frozen=True)
class DirectionalRRL:
    threshold: float
    support_max: float
    status: str
    n_valid_pairs: int


@dataclass(frozen=True)
class PairRRL:
    A_to_B: DirectionalRRL
    B_to_A: DirectionalRRL
    conservative: float
    conservative_status: str
    conservative_lower_bound: float


def directional_rrl(base_scores, perturbed_scores, alpha: float = 0.05) -> DirectionalRRL:
    """Observed-support directional Rank Resolution Limit.

    No extrapolation is performed. If the fitted reversal probability does not
    reach ``alpha`` on sorted unique observed score-gap support, the threshold is
    returned as NaN with status ``NOT_REACHED`` and support_max as a right-censoring
    bound.
    """
    base = np.asarray(base_scores, float)
    pert = np.asarray(perturbed_scores, float)
    if base.shape != pert.shape or base.ndim != 1:
        raise ValueError("base_scores and perturbed_scores must be equal-length 1D arrays")
    if not (0 < alpha < 1):
        raise ValueError("alpha must lie in (0,1)")
    ii, jj = np.triu_indices(len(base), 1)
    d0 = base[ii] - base[jj]
    d1 = pert[ii] - pert[jj]
    valid = (d0 != 0) & (d1 != 0)
    gap = np.abs(d0[valid])
    rev = (np.sign(d0[valid]) != np.sign(d1[valid])).astype(float)
    if len(gap) == 0:
        return DirectionalRRL(np.nan, np.nan, "NO_VALID_PAIRS", 0)
    model = IsotonicRegression(increasing=False, out_of_bounds="clip")
    model.fit(gap, rev)
    support = np.unique(gap)
    support.sort()
    pred = model.predict(support)
    support_max = float(support[-1])
    q = np.where(pred <= alpha)[0]
    if len(q):
        return DirectionalRRL(float(support[q[0]]), support_max, "FINITE", len(gap))
    return DirectionalRRL(np.nan, support_max, "NOT_REACHED", len(gap))


def conservative_pair_rrl(scores_A, scores_B, alpha: float = 0.05) -> PairRRL:
    a = directional_rrl(scores_A, scores_B, alpha=alpha)
    b = directional_rrl(scores_B, scores_A, alpha=alpha)
    if np.isfinite(a.threshold) and np.isfinite(b.threshold):
        cons = max(a.threshold, b.threshold)
        return PairRRL(a, b, cons, "FINITE", cons)
    bounds = [a.threshold if np.isfinite(a.threshold) else a.support_max,
              b.threshold if np.isfinite(b.threshold) else b.support_max]
    return PairRRL(a, b, np.nan, "CENSORED", float(np.nanmax(bounds)))
