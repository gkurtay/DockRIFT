from __future__ import annotations

from itertools import combinations
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, kendalltau

from .ranks import add_seed_ranks, consensus_table
from .rrl import conservative_pair_rrl


def seed_reproducibility(x: pd.DataFrame, receptor: str) -> dict:
    y = add_seed_ranks(x[x.receptor == receptor])
    w = y.pivot(index="ligand_id", columns="seed", values="rank")
    vals = []
    for s1, s2 in combinations(sorted(w.columns), 2):
        vals.append((
            float(spearmanr(w[s1], w[s2]).statistic),
            float(kendalltau(w[s1], w[s2]).statistic),
        ))
    if not vals:
        raise ValueError("At least two seeds are required")
    return {
        "mean_seed_rho": float(np.mean([v[0] for v in vals])),
        "mean_seed_tau": float(np.mean([v[1] for v in vals])),
        "n_seed_pairs": len(vals),
    }


def analyze_receptor_pair(x: pd.DataFrame, A: str, B: str, alpha: float = 0.05,
                          topk=(4,20,40,80)) -> dict:
    A, B = A.lower(), B.lower()
    if A == B:
        raise ValueError("A and B must differ")
    sub = x[x.receptor.isin([A, B])].copy()
    cons = consensus_table(sub)
    smed = cons.pivot(index="ligand_id", columns="receptor", values="score_median")
    ranks = cons.pivot(index="ligand_id", columns="receptor", values="consensus_rank")
    rmin = cons.pivot(index="ligand_id", columns="receptor", values="rank_min")
    rmax = cons.pivot(index="ligand_id", columns="receptor", values="rank_max")
    if A not in smed or B not in smed:
        raise ValueError("Requested receptor absent")
    order = sorted(set(smed.index))
    smed = smed.reindex(order); ranks = ranks.reindex(order); rmin = rmin.reindex(order); rmax = rmax.reindex(order)
    sA, sB = smed[A].to_numpy(float), smed[B].to_numpy(float)
    rA, rB = ranks[A].to_numpy(float), ranks[B].to_numpy(float)
    rho = float(spearmanr(rA, rB).statistic)
    tau = float(kendalltau(rA, rB).statistic)
    abs_drank = np.abs(rA-rB)
    ii, jj = np.triu_indices(len(order), 1)
    dA, dB = sA[ii]-sA[jj], sB[ii]-sB[jj]
    valid = (dA != 0) & (dB != 0)
    inv = float((np.sign(dA[valid]) != np.sign(dB[valid])).mean())
    max_seed_width = np.maximum((rmax[A]-rmin[A]).to_numpy(float),
                                (rmax[B]-rmin[B]).to_numpy(float))
    rec_gt_seed = float(np.mean(abs_drank > max_seed_width))
    out = {
        "pdb_A": A,
        "pdb_B": B,
        "n_ligands": len(order),
        "receptor_spearman_rho": rho,
        "receptor_kendall_tau": tau,
        "median_abs_delta_rank": float(np.median(abs_drank)),
        "mean_abs_delta_rank": float(np.mean(abs_drank)),
        "p90_abs_delta_rank": float(np.quantile(abs_drank, .90)),
        "pair_reversal_fraction": inv,
        "receptor_gt_seed_fraction": rec_gt_seed,
    }
    for k in topk:
        if k > len(order):
            continue
        SA, SB = set(np.argsort(sA)[:k]), set(np.argsort(sB)[:k])
        shared = len(SA & SB)
        out[f"top{k}_survival"] = shared/k
        out[f"top{k}_jaccard"] = shared/len(SA|SB)
    sa, sb = seed_reproducibility(sub, A), seed_reproducibility(sub, B)
    out["mean_seed_rho_A"] = sa["mean_seed_rho"]
    out["mean_seed_rho_B"] = sb["mean_seed_rho"]
    out["mean_pair_seed_rho"] = float(np.mean([sa["mean_seed_rho"], sb["mean_seed_rho"]]))
    out["mean_seed_tau_A"] = sa["mean_seed_tau"]
    out["mean_seed_tau_B"] = sb["mean_seed_tau"]
    out["mean_pair_seed_tau"] = float(np.mean([sa["mean_seed_tau"], sb["mean_seed_tau"]]))
    rrl = conservative_pair_rrl(sA, sB, alpha=alpha)
    out.update({
        "rrl_alpha": alpha,
        "rrl_A_to_B": rrl.A_to_B.threshold,
        "rrl_A_status": rrl.A_to_B.status,
        "rrl_A_support_max": rrl.A_to_B.support_max,
        "rrl_B_to_A": rrl.B_to_A.threshold,
        "rrl_B_status": rrl.B_to_A.status,
        "rrl_B_support_max": rrl.B_to_A.support_max,
        "rrl_conservative": rrl.conservative,
        "rrl_conservative_status": rrl.conservative_status,
        "rrl_conservative_lower_bound": rrl.conservative_lower_bound,
    })
    return out
