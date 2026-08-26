from pathlib import Path
import pandas as pd
from dockrift.io import read_score_table
from dockrift.report import analyze_score_file
from dockrift.cli import doctor
from dockrift.gui.logic import load_demo

def test_multiscoring_cli_reader(tmp_path):
    s,_,_=load_demo();p=tmp_path/"scores.csv";s.data.to_csv(p,index=False)
    x=read_score_table(p);assert set(x.scoring_function)=={"vina","vinardo"}
    out=analyze_score_file(p);assert len(out)==2

def test_doctor_ready():
    assert doctor()["ready"] is True
