import math
from dockrift.gui.logic import load_demo, pair_metrics

def test_vdr_regression_vina():
    s,_,_=load_demo();m=pair_metrics(s.data,"vina","3a2i","3a2j",.05,20)
    assert math.isclose(m["receptor_spearman_rho"],0.9585856068319556,abs_tol=1e-12)
    assert math.isclose(m["mean_pair_seed_rho"],0.9974781596445379,abs_tol=1e-12)
    assert m["rrl_display"]["status"]=="FINITE"
    assert math.isclose(m["rrl_display"]["value"],0.875,abs_tol=1e-12)

def test_vdr_regression_vinardo():
    s,_,_=load_demo();m=pair_metrics(s.data,"vinardo","3a2i","3a2j",.05,20)
    assert math.isclose(m["receptor_spearman_rho"],0.9517698268798475,abs_tol=1e-12)
    assert math.isclose(m["rrl_display"]["value"],0.8840000000000012,abs_tol=1e-12)
