from __future__ import annotations
import numpy as np
import pandas as pd


def normalize_scores(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize common docking-score table aliases to DockRIFT long-form schema."""
    x=df.copy(); lower={str(c).lower():c for c in x.columns}
    aliases={
        'ligand_id':['ligand_id','dockrift_id','ligand','compound_id','id'],
        'receptor':['receptor','pdb','pdb_id','structure','conformation'],
        'seed':['seed','random_seed','run_seed'],
        'score':['score','mode1_score','vina_score','affinity'],
        'scoring_function':['scoring_function','scoring','engine'],
        'class':['class','activity_class','label','activity'],
    }
    rename={}
    for canonical,opts in aliases.items():
        for opt in opts:
            if opt in lower:
                rename[lower[opt]]=canonical; break
    x=x.rename(columns=rename)
    required=['ligand_id','receptor','seed','score']
    missing=[c for c in required if c not in x.columns]
    if missing: raise ValueError('Missing required columns after automatic mapping: '+', '.join(missing))
    if 'scoring_function' not in x.columns: x['scoring_function']='vina'
    keep=required+['scoring_function']+(['class'] if 'class' in x.columns else [])
    x=x[keep].copy(); x['ligand_id']=x['ligand_id'].astype(str)
    x['receptor']=x['receptor'].astype(str).str.strip().str.lower(); x['scoring_function']=x['scoring_function'].astype(str).str.strip().str.lower()
    x['seed']=pd.to_numeric(x['seed'],errors='raise').astype(int); x['score']=pd.to_numeric(x['score'],errors='raise').astype(float)
    if 'class' in x.columns:
        x['class']=x['class'].astype(str).str.strip().str.lower().replace({'1':'active','true':'active','yes':'active','positive':'active','0':'inactive','false':'inactive','no':'inactive','negative':'inactive'})
    keys=['scoring_function','receptor','seed','ligand_id']
    if x.duplicated(keys).any(): raise ValueError('Duplicate scoring/receptor/seed/ligand keys detected')
    if not np.isfinite(x['score']).all(): raise ValueError('Non-finite score values detected')
    return x.sort_values(keys,kind='stable').reset_index(drop=True)
