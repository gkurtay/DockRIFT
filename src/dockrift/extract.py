from __future__ import annotations
from pathlib import Path
import re
import pandas as pd

REMARK_RE=re.compile(r"^REMARK\s+VINA\s+RESULT:\s+([-+0-9.eE]+)",re.I)


def mode1_score_from_pdbqt(path:str|Path)->float:
    p=Path(path)
    if not p.is_file() or p.stat().st_size==0:raise FileNotFoundError(f"Missing/empty PDBQT: {p}")
    for line in p.read_text(errors="replace").splitlines():
        m=REMARK_RE.match(line.strip())
        if m:return float(m.group(1))
    raise ValueError(f"No REMARK VINA RESULT score found: {p}")


def score_table_from_manifest(manifest_path:str|Path,out_csv:str|Path|None=None)->pd.DataFrame:
    p=Path(manifest_path);m=pd.read_csv(p,keep_default_na=False)
    required={"scoring_function","receptor","ligand_id","seed","out_pdbqt"};missing=required-set(m.columns)
    if missing:raise ValueError(f"Manifest missing required columns: {sorted(missing)}")
    key=["scoring_function","receptor","ligand_id","seed"]
    if m[key].duplicated().any():raise ValueError("Manifest contains duplicate scoring/receptor/ligand/seed keys")
    rows=[]
    for _,r in m.iterrows():
        row={"scoring_function":str(r["scoring_function"]).lower(),"receptor":str(r["receptor"]).lower(),"ligand_id":str(r["ligand_id"]),"seed":int(r["seed"]),"score":mode1_score_from_pdbqt(r["out_pdbqt"])}
        if "class" in m.columns:row["class"]=str(r["class"])
        rows.append(row)
    x=pd.DataFrame(rows).sort_values(["scoring_function","receptor","seed","ligand_id"],kind="stable").reset_index(drop=True)
    if out_csv is not None:
        q=Path(out_csv);q.parent.mkdir(parents=True,exist_ok=True);x.to_csv(q,index=False)
    return x
