from pathlib import Path
import zipfile
from dockrift.gui.logic import load_demo,dataset_profile,data_health,pair_payload,seed_convergence,validation_panel,export_bundle,provenance_payload,bootstrap_panel

def test_profile_and_health():
    s,_,_=load_demo();p=dataset_profile(s.data,s.source_name,s.sha256);h=data_health(s.data)
    assert p["rows"]==4800 and p["ligands"]==400 and p["n_receptors"]==2 and p["n_seeds"]==3
    assert p["balanced"] and p["has_labels"]
    assert len(h["table"])==12

def test_pair_payload_has_all_studio_views():
    s,_,_=load_demo();p=pair_payload(s.data,"vina","3a2i","3a2j",.05,20)
    for key in ["metrics","rank_figure","delta_figure","rrl_figure","alpha_figure","topk_figure","cross_scoring"]: assert key in p
    assert len(p["cross_scoring"])==2

def test_seed_convergence_and_validation():
    s,_,_=load_demo();c=seed_convergence(s.data,"vina","3a2i","3a2j");v=validation_panel(s.data,"vina","3a2i","3a2j",20)
    assert [r["n_seeds"] for r in c["table"]]==[2,3]
    assert "abs_delta_roc_auc" in v["metrics"]

def test_export_bundle():
    s,_,_=load_demo();p=Path(export_bundle(s,"vina","3a2i","3a2j",.05,20));assert p.is_file()
    with zipfile.ZipFile(p) as z:
        names=set(z.namelist())
    assert {"normalized_scores.csv","pair_metrics.csv","methods_capsule.txt","provenance.json","interactive_report.html"}.issubset(names)

def test_provenance():
    s,_,_=load_demo();p=provenance_payload(s);assert p["source_sha256"]==s.sha256 and p["method"]["rrl"].startswith("directional")


def test_bootstrap_panel():
    s,_,_=load_demo();b=bootstrap_panel(s.data,"vina","3a2i","3a2j",.05,20,20260813)
    assert b["summary"]["bootstrap_finite_n"] + b["summary"]["bootstrap_censored_n"] == 20
    assert "figure" in b
