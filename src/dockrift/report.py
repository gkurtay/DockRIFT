from __future__ import annotations
from itertools import combinations
from pathlib import Path
import pandas as pd
from .io import read_score_table, validate_balanced_design
from .metrics import analyze_receptor_pair


def analyze_score_file(path,*,alpha=.05,scoring_function=None):
    x=read_score_table(path);rows=[]
    scoring_values=[str(scoring_function).lower()] if scoring_function is not None else sorted(x["scoring_function"].unique())
    for scoring in scoring_values:
        q=x[x["scoring_function"]==scoring].copy()
        if q.empty:raise ValueError(f"Scoring function {scoring!r} not present")
        receptors,seeds,n_ligands=validate_balanced_design(q)
        for A,B in combinations(receptors,2):
            r=analyze_receptor_pair(q,A,B,alpha=alpha);r["scoring_function"]=scoring;r["n_seeds"]=len(seeds);rows.append(r)
    return pd.DataFrame(rows)


def write_analysis(path,out_csv,*,alpha=.05,scoring_function=None):
    df=analyze_score_file(path,alpha=alpha,scoring_function=scoring_function);Path(out_csv).parent.mkdir(parents=True,exist_ok=True);df.to_csv(out_csv,index=False);return df
