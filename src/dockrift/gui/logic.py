from __future__ import annotations
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import base64, hashlib, io, json, math, tempfile, zipfile, zlib
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from .. import __version__
from ..bootstrap import bootstrap_pair_rrl, summarize_bootstrap
from ..metrics import analyze_receptor_pair
from ..ranks import consensus_table
from ..validation import pair_validity_instability, screening_metrics, heldout_rrl_calibration, subsample_rrl_sensitivity
from ..studio_utils import normalize_scores

@dataclass
class SessionData:
    data: pd.DataFrame
    source_name: str
    sha256: str
    source_kind: str="upload"

def _sha(b:bytes)->str:return hashlib.sha256(b).hexdigest()
def _json_ready(x):
    if isinstance(x,dict):return {k:_json_ready(v) for k,v in x.items()}
    if isinstance(x,list):return [_json_ready(v) for v in x]
    if isinstance(x,(np.integer,)):return int(x)
    if isinstance(x,(np.floating,)):x=float(x)
    if isinstance(x,float) and not math.isfinite(x):return None
    return x

def _fig_json(fig):return json.loads(fig.to_json())
def read_score_bytes(filename:str,content:bytes)->SessionData:
    sep='\t' if str(filename).lower().endswith(('.tsv','.txt')) else ','
    raw=pd.read_csv(io.BytesIO(content),sep=sep);x=normalize_scores(raw)
    return SessionData(x,Path(filename).name,_sha(content),'upload')

def load_demo():
    demo=Path(__file__).resolve().parents[1]/'demo'
    source=demo/'vdr_vina_vinardo_top1_scores.csv'
    if source.is_file():
        raw=pd.read_csv(source)
        digest_bytes=source.read_bytes()
    else:
        numbered=[]
        for q in demo.glob('vdr_demo_b85_part*.txt'):
            try:numbered.append((int(q.stem.rsplit('part',1)[1]),q))
            except Exception:pass
        parts=[q for _,q in sorted(numbered)]
        if not parts:raise FileNotFoundError('Packaged VDR demo table is missing')
        payload=''.join(q.read_text(encoding='ascii').strip() for q in parts).encode('ascii')
        digest_bytes=zlib.decompress(base64.b85decode(payload))
        raw=pd.read_csv(io.BytesIO(digest_bytes))
    x=normalize_scores(raw)
    master=pd.read_csv(demo/'dockrift_publication_pair_metric_master.csv')
    bio=pd.read_csv(demo/'dockrift_final_biological_target_scoring_summary.csv')
    return SessionData(x,'Built-in VDR dual-scoring demo',_sha(digest_bytes),'demo'),master,bio

def dataset_profile(x,source_name,sha256):
    sc=sorted(x.scoring_function.unique().tolist());rec=sorted(x.receptor.unique().tolist());seeds=sorted(int(v) for v in x.seed.unique());n=int(x.ligand_id.nunique())
    per=x.groupby(['scoring_function','receptor','seed']).ligand_id.nunique();balanced=bool(len(per) and per.nunique()==1 and int(per.iloc[0])==n)
    labels='class' in x.columns and bool(set(x['class'].unique())&{'active','inactive'})
    return {'source_name':source_name,'sha256':sha256,'rows':int(len(x)),'ligands':n,'receptors':rec,'n_receptors':len(rec),'seeds':seeds,'n_seeds':len(seeds),'scoring_functions':sc,'receptors_by_scoring':{s:sorted(x.loc[x.scoring_function==s,'receptor'].unique().tolist()) for s in sc},'balanced':balanced,'has_labels':labels,'ligand_ids':sorted(x.ligand_id.unique().tolist())}

def data_health(x):
    q=x.groupby(['scoring_function','receptor','seed'],as_index=False).agg(rows=('ligand_id','size'),ligands=('ligand_id','nunique'),score_min=('score','min'),score_median=('score','median'),score_max=('score','max'))
    fig=go.Figure(go.Bar(x=[f"{r.scoring_function}:{r.receptor}:{r.seed}" for _,r in q.iterrows()],y=q.ligands));fig.update_layout(title='Design coverage',height=360)
    return {'table':q.to_dict(orient='records'),'coverage_figure':_fig_json(fig)}

def _consensus(x,sc,A,B):
    sub=x[(x.scoring_function==sc.lower())&x.receptor.isin([A.lower(),B.lower()])].copy();c=consensus_table(sub);return sub,c

def _rrl_display(m):
    if m.get('rrl_conservative_status')=='FINITE' and np.isfinite(m.get('rrl_conservative',np.nan)):return {'status':'FINITE','value':float(m['rrl_conservative']),'label':f"{m['rrl_conservative']:.3f}"}
    lb=m.get('rrl_conservative_lower_bound',np.nan);return {'status':'RIGHT_CENSORED','value':float(lb) if np.isfinite(lb) else None,'label':f"≥ {lb:.3f}" if np.isfinite(lb) else 'not estimable'}

def pair_metrics(x,sc,A,B,alpha=.05,topk=20):
    m=analyze_receptor_pair(x[x.scoring_function==sc.lower()].copy(),A,B,alpha=float(alpha),topk=tuple(sorted(set([4,int(topk),40,80]))));m['transfer_gap']=float(m['mean_pair_seed_rho']-m['receptor_spearman_rho']);m['rrl_display']=_rrl_display(m);return _json_ready(m)

def _rank_fig(x,sc,A,B):
    _,c=_consensus(x,sc,A,B);p=c.pivot(index='ligand_id',columns='receptor',values='consensus_rank');fig=go.Figure(go.Scattergl(x=p[A.lower()],y=p[B.lower()],mode='markers'));n=len(p);fig.add_shape(type='line',x0=1,y0=1,x1=n,y1=n,line={'dash':'dash'});fig.update_layout(title='Rank-transfer map',xaxis_title=A.upper(),yaxis_title=B.upper(),height=470);return fig

def _delta_fig(x,sc,A,B):
    _,c=_consensus(x,sc,A,B);p=c.pivot(index='ligand_id',columns='receptor',values='consensus_rank');d=(p[B.lower()]-p[A.lower()]).sort_values(key=lambda s:s.abs(),ascending=False).head(30);fig=go.Figure(go.Bar(x=d.values,y=d.index,orientation='h'));fig.update_layout(title='Largest rank displacements',height=560);return fig

def _rrl_fig(x,sc,A,B,alpha):
    m=pair_metrics(x,sc,A,B,alpha,20);vals=[m.get('rrl_A_to_B'),m.get('rrl_B_to_A')];fig=go.Figure(go.Bar(x=[f'{A}→{B}',f'{B}→{A}'],y=[v or 0 for v in vals]));fig.update_layout(title='Directional RRL',height=360);return fig

def pair_payload(x,sc,A,B,alpha=.05,topk=20):
    m=pair_metrics(x,sc,A,B,alpha,topk);cross=[]
    for s in sorted(x.scoring_function.unique()):
        if {A.lower(),B.lower()}.issubset(set(x.loc[x.scoring_function==s,'receptor'].unique())):cross.append({'scoring_function':s,**pair_metrics(x,s,A,B,alpha,topk)})
    vals=[]
    for a in [.01,.025,.05,.1]:vals.append(pair_metrics(x,sc,A,B,a,topk)['rrl_display']['value'])
    alpha_fig=go.Figure(go.Scatter(x=[.01,.025,.05,.1],y=vals,mode='lines+markers'));alpha_fig.update_layout(title='RRL criterion sensitivity')
    topfig=go.Figure(go.Bar(x=['Top-k survival'],y=[m.get(f'top{int(topk)}_survival',0)]));topfig.update_layout(title='Shortlist stability')
    return {'metrics':m,'rank_figure':_fig_json(_rank_fig(x,sc,A,B)),'delta_figure':_fig_json(_delta_fig(x,sc,A,B)),'rrl_figure':_fig_json(_rrl_fig(x,sc,A,B,alpha)),'alpha_figure':_fig_json(alpha_fig),'topk_figure':_fig_json(topfig),'cross_scoring':cross}

def seed_convergence(x,sc,A,B):
    q=x[x.scoring_function==sc.lower()].copy();seeds=sorted(q.seed.unique());rows=[]
    for n in range(2,len(seeds)+1):
        z=q[q.seed.isin(seeds[:n])];m=analyze_receptor_pair(z,A,B);rows.append({'n_seeds':n,'mean_pair_seed_rho':m['mean_pair_seed_rho'],'receptor_spearman_rho':m['receptor_spearman_rho']})
    fig=go.Figure();fig.add_trace(go.Scatter(x=[r['n_seeds'] for r in rows],y=[r['mean_pair_seed_rho'] for r in rows],name='seed ρ'));fig.add_trace(go.Scatter(x=[r['n_seeds'] for r in rows],y=[r['receptor_spearman_rho'] for r in rows],name='receptor ρ'));return {'table':rows,'figure':_fig_json(fig)}

def validation_panel(x,sc,A,B,topk=20):
    sub,c=_consensus(x,sc,A,B)
    if 'class' not in sub.columns:return {'metrics':{},'figure':_fig_json(go.Figure())}
    p=c.pivot(index='ligand_id',columns='receptor',values='score_median');lab=sub.groupby('ligand_id')['class'].first().reindex(p.index);y=(lab=='active').astype(int).to_numpy();m=pair_validity_instability(y,p[A.lower()].to_numpy(),p[B.lower()].to_numpy(),topk=int(topk));fig=go.Figure(go.Bar(x=['ROC-AUC Δ','PR-AUC Δ',f'EF@{topk} Δ'],y=[m['abs_delta_roc_auc'],m['abs_delta_pr_auc'],m['abs_delta_ef_topk']]));return {'metrics':_json_ready(m),'figure':_fig_json(fig)}

def provenance_payload(s):return {'dockrift_version':__version__,'source_name':s.source_name,'source_sha256':s.sha256,'source_kind':s.source_kind,'method':{'rrl':'directional observed-support isotonic reversal threshold; no extrapolation; explicit censoring','ranking':'ascending score; median-across-seed consensus'}}

def export_bundle(s,sc,A,B,alpha=.05,topk=20):
    tmp=Path(tempfile.mkdtemp(prefix='dockrift_export_'));zpath=tmp/'dockrift_reproducibility_export.zip';m=pair_metrics(s.data,sc,A,B,alpha,topk);report=f"<html><body><h1>DockRIFT {__version__}</h1><p>{A.upper()} ↔ {B.upper()} · {sc}</p><pre>{json.dumps(m,indent=2)}</pre></body></html>"
    with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED) as z:
        z.writestr('normalized_scores.csv',s.data.to_csv(index=False));z.writestr('pair_metrics.csv',pd.DataFrame([m]).to_csv(index=False));z.writestr('methods_capsule.txt','DockRIFT observed-support RRL; no extrapolation; censoring preserved.\n');z.writestr('provenance.json',json.dumps(provenance_payload(s),indent=2));z.writestr('interactive_report.html',report)
    return str(zpath)

def bootstrap_panel(x,sc,A,B,alpha=.05,n_boot=20,rng_seed=20260813):
    _,c=_consensus(x,sc,A,B);p=c.pivot(index='ligand_id',columns='receptor',values='score_median');rows=bootstrap_pair_rrl(p[A.lower()].to_numpy(),p[B.lower()].to_numpy(),n_boot=int(n_boot),rng_seed=int(rng_seed),alpha=float(alpha));s=summarize_bootstrap(rows,int(n_boot));fig=go.Figure(go.Histogram(x=[r['rrl_conservative'] for r in rows if np.isfinite(r['rrl_conservative'])]));return {'summary':_json_ready(s),'table':_json_ready(rows),'figure':_fig_json(fig)}
def robustness_panel(x,sc,A,B,alpha=.05,heldout_repeats=30,subsample_repeats=20):
    _,c=_consensus(x,sc,A,B);p=c.pivot(index='ligand_id',columns='receptor',values='score_median');a=p[A.lower()].to_numpy();b=p[B.lower()].to_numpy();h=heldout_rrl_calibration(a,b,repeats=int(heldout_repeats),alpha=float(alpha));sizes=[n for n in (100,200,300) if n<=len(a)];u=subsample_rrl_sensitivity(a,b,sample_sizes=sizes,repeats=int(subsample_repeats),alpha=float(alpha));return {'heldout':_json_ready(h),'universe':_json_ready(u)}
def pair_matrix_figure(x,sc,metric='receptor_spearman_rho'):
    rec=sorted(x.loc[x.scoring_function==sc.lower(),'receptor'].unique());rows=[]
    for A,B in combinations(rec,2):rows.append(pair_metrics(x,sc,A,B))
    fig=go.Figure(go.Bar(x=[f"{r['pdb_A']}-{r['pdb_B']}" for r in rows],y=[r.get(metric,0) for r in rows]));return fig,rows
def seed_vs_receptor_figure(x,sc):
    rec=sorted(x.loc[x.scoring_function==sc.lower(),'receptor'].unique());rows=[pair_metrics(x,sc,A,B) for A,B in combinations(rec,2)];fig=go.Figure(go.Scatter(x=[r['mean_pair_seed_rho'] for r in rows],y=[r['receptor_spearman_rho'] for r in rows],mode='markers'));return fig
def screening_landscape(x,sc,topk=20):
    rows=[]
    if 'class' not in x.columns:return {'table':rows}
    for r in sorted(x.loc[x.scoring_function==sc.lower(),'receptor'].unique()):
        q=x[(x.scoring_function==sc.lower())&(x.receptor==r)];c=consensus_table(q);p=c.set_index('ligand_id')['score_median'];lab=q.groupby('ligand_id')['class'].first().reindex(p.index);rows.append({'receptor':r,**screening_metrics((lab=='active').astype(int).to_numpy(),p.to_numpy(),topk=int(topk))})
    return {'table':_json_ready(rows)}
def ligand_inspector(x,sc,ligand_id):
    q=x[(x.scoring_function==sc.lower())&(x.ligand_id.astype(str)==str(ligand_id))].sort_values(['receptor','seed']);fig=go.Figure(go.Scatter(x=[f"{r.receptor}:{r.seed}" for _,r in q.iterrows()],y=q.score,mode='lines+markers'));return {'table':q.to_dict(orient='records'),'figure':_fig_json(fig)}
def demo_atlas(master,bio):return {'master_rows':int(len(master)),'biological_rows':int(len(bio)),'master':master.to_dict(orient='records'),'biological':bio.to_dict(orient='records')}
