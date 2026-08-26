from __future__ import annotations

import numpy as np
from sklearn.isotonic import IsotonicRegression


def _direction_from_pair_arrays(gap, rev, alpha):
    if len(gap) == 0:
        return np.nan, np.nan
    m = IsotonicRegression(increasing=False, out_of_bounds="clip")
    m.fit(gap, rev)
    support = np.unique(gap); support.sort()
    pred = m.predict(support)
    mx = float(support[-1])
    q = np.where(pred <= alpha)[0]
    return (float(support[q[0]]), mx) if len(q) else (np.nan, mx)


def bootstrap_pair_rrl(scores_A, scores_B, *, n_boot=300, rng_seed=20260813,
                       alpha=0.05):
    """Ligand-identity bootstrap preserving the frozen DockRIFT duplicate rule."""
    A = np.asarray(scores_A, float); B = np.asarray(scores_B, float)
    if len(A) != len(B):
        raise ValueError("Score arrays differ in length")
    n = len(A)
    tri_i, tri_j = np.triu_indices(n, 1)
    rng = np.random.default_rng(rng_seed)
    rows = []
    for b in range(n_boot):
        sampled = rng.integers(0, n, size=n)
        a, c = A[sampled], B[sampled]
        distinct = sampled[tri_i] != sampled[tri_j]
        dA, dB = a[tri_i]-a[tri_j], c[tri_i]-c[tri_j]
        valid = distinct & (dA != 0) & (dB != 0)
        rev = (np.sign(dA[valid]) != np.sign(dB[valid])).astype(float)
        ta, maxa = _direction_from_pair_arrays(np.abs(dA[valid]), rev, alpha)
        tb, maxb = _direction_from_pair_arrays(np.abs(dB[valid]), rev, alpha)
        if np.isfinite(ta) and np.isfinite(tb):
            cons, status = max(ta, tb), "FINITE"
        else:
            cons, status = np.nan, "CENSORED"
        rows.append({
            "bootstrap": b+1,
            "rrl_A_to_B": ta,
            "rrl_A_support_max": maxa,
            "rrl_B_to_A": tb,
            "rrl_B_support_max": maxb,
            "rrl_conservative": cons,
            "status": status,
        })
    return rows


def summarize_bootstrap(rows, n_boot=None):
    vals = np.asarray([r["rrl_conservative"] for r in rows
                       if np.isfinite(r["rrl_conservative"])], float)
    n_boot = len(rows) if n_boot is None else n_boot
    censored = n_boot - len(vals)
    if len(vals) == n_boot:
        status = "FULLY_ESTIMABLE"
    elif len(vals):
        status = "FINITE_ONLY_DIAGNOSTIC"
    else:
        status = "NO_FINITE_BOOTSTRAPS"
    if len(vals):
        lo, med, hi = np.percentile(vals, [2.5, 50, 97.5])
    else:
        lo = med = hi = np.nan
    return {
        "bootstrap_finite_n": len(vals),
        "bootstrap_censored_n": censored,
        "bootstrap_finite_only_median": med,
        "bootstrap_finite_only_lo": lo,
        "bootstrap_finite_only_hi": hi,
        "bootstrap_interval_status": status,
    }
