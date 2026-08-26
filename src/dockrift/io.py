from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

ALIASES={
    "dockrift_id":"ligand_id","dock_id":"ligand_id","ligand":"ligand_id","compound_id":"ligand_id","id":"ligand_id",
    "pdb_id":"receptor","pdb":"receptor","structure":"receptor","conformation":"receptor",
    "mode1_score":"score","score_kcal_mol":"score","vina_score":"score","affinity":"score",
    "scoring":"scoring_function","engine":"scoring_function","activity_class":"class","label":"class","activity":"class",
}
REQUIRED={"ligand_id","receptor","seed","score"}


def read_score_table(path: str|Path)->pd.DataFrame:
    p=Path(path);sep="\t" if p.suffix.lower() in {".tsv",".txt"} else ",";x=pd.read_csv(p,sep=sep,keep_default_na=False)
    lower={c.lower():c for c in x.columns};rename={}
    for src,dst in ALIASES.items():
        if dst not in x.columns and src in lower: rename[lower[src]]=dst
    x=x.rename(columns=rename);missing=REQUIRED-set(x.columns)
    if missing: raise ValueError(f"Missing required score-table columns: {sorted(missing)}")
    if "scoring_function" not in x.columns:x["scoring_function"]="vina"
    keep=["ligand_id","receptor","seed","score","scoring_function"]+(["class"] if "class" in x.columns else [])
    x=x[keep].copy();x["ligand_id"]=x["ligand_id"].astype(str);x["receptor"]=x["receptor"].astype(str).str.strip().str.lower();x["scoring_function"]=x["scoring_function"].astype(str).str.strip().str.lower();x["seed"]=pd.to_numeric(x["seed"],errors="raise").astype(int);x["score"]=pd.to_numeric(x["score"],errors="raise").astype(float)
    if "class" in x.columns:
        x["class"]=x["class"].astype(str).str.strip().str.lower().replace({"1":"active","true":"active","yes":"active","positive":"active","0":"inactive","false":"inactive","no":"inactive","negative":"inactive"})
    key=["scoring_function","ligand_id","receptor","seed"]
    if x[key].duplicated().any():raise ValueError("Duplicate scoring/ligand/receptor/seed keys in score table")
    if not np.isfinite(x["score"]).all():raise ValueError("Non-finite scores detected")
    return x.sort_values(["scoring_function","receptor","seed","ligand_id"],kind="stable").reset_index(drop=True)


def validate_balanced_design(x:pd.DataFrame)->tuple[list[str],list[int],int]:
    receptors=sorted(x["receptor"].unique().tolist());seeds=sorted(int(v) for v in x["seed"].unique());ligands=sorted(x["ligand_id"].unique().tolist())
    expected=len(receptors)*len(seeds)*len(ligands)
    if len(x)!=expected:raise ValueError(f"Unbalanced design: rows={len(x)}, expected={expected} for {len(receptors)} receptors x {len(seeds)} seeds x {len(ligands)} ligands")
    for rec in receptors:
        for seed in seeds:
            y=x[(x.receptor==rec)&(x.seed==seed)]
            if set(y.ligand_id)!=set(ligands):raise ValueError(f"Ligand identities differ for receptor={rec}, seed={seed}")
    return receptors,seeds,len(ligands)
