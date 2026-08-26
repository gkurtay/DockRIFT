from pathlib import Path
import pandas as pd
from dockrift.extract import mode1_score_from_pdbqt,score_table_from_manifest
from dockrift.provenance import write_sha256_manifest,verify_sha256_manifest

def test_extract_manifest(tmp_path):
    out=tmp_path/"x.pdbqt";out.write_text("MODEL 1\nREMARK VINA RESULT: -8.123 0.0 0.0\nENDMDL\n")
    m=pd.DataFrame([{"scoring_function":"vina","receptor":"1abc","ligand_id":"L1","class":"active","seed":101,"out_pdbqt":str(out)}]);mp=tmp_path/"manifest.csv";m.to_csv(mp,index=False)
    x=score_table_from_manifest(mp,tmp_path/"scores.csv");assert len(x)==1 and x.iloc[0].score==-8.123 and x.iloc[0]["class"]=="active"

def test_sha_freeze(tmp_path):
    a=tmp_path/"a.txt";a.write_text("DockRIFT\n");m=write_sha256_manifest([a],tmp_path/"freeze.sha256");rows=verify_sha256_manifest(m);assert rows[0]["ok"]
    a.write_text("changed\n");assert not verify_sha256_manifest(m)[0]["ok"]
