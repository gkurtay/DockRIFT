from __future__ import annotations
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import hashlib, io, json, math, platform, sys, tempfile, zipfile
from typing import Any
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.stats import spearmanr
from .. import __version__
from ..bootstrap import bootstrap_pair_rrl, summarize_bootstrap
from ..metrics import analyze_receptor_pair, seed_reproducibility
from ..ranks import consensus_table
from ..rrl import conservative_pair_rrl
from ..validation import heldout_rrl_calibration, pair_validity_instability, screening_metrics, subsample_rrl_sensitivity
from ..studio_utils import normalize_scores

BRAND={'cyan':'#43d9ea','blue':'#5c8dff','magenta':'#f071c9','violet':'#a978ff','text':'#edf5ff','muted':'#8fa7c7','grid':'rgba(143,167,199,.13)'}

@dataclass
class SessionData:
    data: pd.DataFrame; source_name: str; sha256: str; source_kind: str='upload'

def _sha256_bytes(b:bytes)->str: return hashlib.sha256(b).hexdigest()
def _fig_json(fig:go.Figure)->dict: return json.loads(fig.to_json())
def _json_ready(o:Any)->Any:
    if isinstance(o,dict): return {str(k):_json_ready(v) for k,v in o.items()}
    if isinstance(o,(list,tuple)): return [_json_ready(v) for v in o]
    if isinstance(o,np.integer): return int(o)
    if isinstance(o,np.floating): return None if not math.isfinite(float(o)) else float(o)
    if isinstance(o,float): return None if not math.isfinite(o) else o
    return o

def _layout(fig,title,height=430):
    fig.update_layout(title=title,height=height,paper_bgcolor='rgba(0,0,0,0)',plot_bgcolor='rgba(6,17,31,.38)',font=dict(color=BRAND['text']),margin=dict(l=50,r=30,t=55,b=45))
    fig.update_xaxes(gridcolor=BRAND['grid']); fig.update_yaxes(gridcolor=BRAND['grid']); return fig

def read_score_bytes(filename:str,content:bytes)->SessionData:
    sep='\t' if str(filename).lower().endswith(('.tsv','.txt')) else ','
    x=normalize_scores(pd.read_csv(io.BytesIO(content),sep=sep)); return SessionData(x,Path(filename).name,_sha256_bytes(content),'upload')

def load_demo():
    demo=Path(__file__).resolve().parents[1]/'demo'; files=sorted((demo/'vdr_branches').glob('*.csv'))
    if len(files)!=12: raise FileNotFoundError('Packaged VDR demo branches are missing')
    raw=pd.concat([pd.read_csv(p) for p in files],ignore_index=True); x=normalize_scores(raw)
    digest=raw.sort_values(['scoring_function','receptor','seed','ligand_id'],kind='stable').to_csv(index=False).encode()
    return SessionData(x,'Built-in VDR dual-scoring demo',_sha256_bytes(digest),'demo'),pd.read_csv(demo/'dockrift_publication_pair_metric_master.csv'),pd.read_csv(demo/'dockrift_final_biological_target_scoring_summary.csv')

def dataset_profile(x,source_name,sha256):
    recs=sorted(x.receptor.unique()); seeds=sorted(map(int,x.seed.unique())); scores=sorted(x.scoring_function.unique()); n=int(x.ligand_id.nunique())
    per=x.groupby(['scoring_function','receptor','seed']).ligand_id.nunique(); balanced=bool(len(per) and per.nunique()==1 and int(per.iloc[0])==n)
    has='class' in x and bool(set(x['class'].dropna()).intersection({'active','inactive'}))
    return {'source_name':source_name,'sha256':sha256,'rows':int(len(x)),'ligands':n,'receptors':recs,'n_receptors':len(recs),'seeds':seeds,'n_seeds':len(seeds),'scoring_functions':scores,'receptors_by_scoring':{s:sorted(x.loc[x.scoring_function==s,'receptor'].unique()) for s in scores},'balanced':balanced,'has_labels':has,'active_ligands':int(x.loc[x.get('class','')=='active','ligand_id'].nunique()) if 'class' in x else 0,'inactive_ligands':int(x.loc[x.get('class','')=='inactive','ligand_id'].nunique()) if 'class' in x else 0,'ligand_ids':sorted(x.ligand_id.astype(str).unique()),'warnings':[],'preview':x.head(20).to_dict('records')}

def data_health(x):
    q=x.groupby(['scoring_function','receptor','seed'],as_index=False).agg(rows=('ligand_id','size'),ligands=('ligand_id','nunique'),score_min=('score','min'),score_median=('score','median'),score_max=('score','max'))
    piv=q.assign(branch=q.scoring_function.str.title()+' · '+q.receptor.str.upper()).pivot(index='branch',columns='seed',values='ligands')
    fig=go.Figure(go.Heatmap(z=piv.to_numpy(float),x=[str(c) for c in piv.columns],y=piv.index.tolist(),colorscale='Viridis'))
    return {'table':q.to_dict('records'),'coverage_figure':_fig_json(_layout(fig,'Design coverage matrix',380))}

def pair_metrics(x,scoring,A,B,alpha=.05,topk=20):
    m=analyze_receptor_pair(x[x.scoring_function==scoring.lower()],A,B,alpha=float(alpha),topk=tuple(sorted(set([4,int(topk),40,80]))))
    m['transfer_gap']=float(m['mean_pair_seed_rho']-m['receptor_spearman_rho']); status=str(m['rrl_conservative_status']).upper()
    if status=='FINITE' and np.isfinite(m.get('rrl_conservative',np.nan)): d={'status':'FINITE','value':float(m['rrl_conservative']),'label':f"{m['rrl_conservative']:.3f}"}
    else:
        lb=m.get('rrl_conservative_lower_bound',np.nan); d={'status':'RIGHT_CENSORED' if np.isfinite(lb) else status,'value':float(lb) if np.isfinite(lb) else None,'label':f'≥ {lb:.3f}' if np.isfinite(lb) else 'not estimable'}
    m['rrl_display']=d; m['topk']=int(topk); return _json_ready(m)

def _pair_consensus(x,scoring,A,B):
    sub=x[(x.scoring_function==scoring.lower()) & x.receptor.isin([A.lower(),B.lower()])]; cons=consensus_table(sub)
    return sub,cons.pivot(index='ligand_id',columns='receptor',values='consensus_rank'),cons.pivot(index='ligand_id',columns='receptor',values='score_median')

def rank_transfer_figure(x,scoring,A,B):
    _,r,_=_pair_consensus(x,scoring,A,B); a,b=A.lower(),B.lower(); fig=go.Figure(go.Scattergl(x=r[a],y=r[b],mode='markers',marker=dict(size=5,color=BRAND['cyan'],opacity=.7)))
    n=len(r); fig.add_trace(go.Scatter(x=[1,n],y=[1,n],mode='lines',line=dict(dash='dash'),showlegend=False)); return _layout(fig,'Rank-transfer map',480)

def delta_rank_figure(x,scoring,A,B):
    _,r,_=_pair_consensus(x,scoring,A,B); d=(r[B.lower()]-r[A.lower()]).sort_values(key=lambda s:s.abs(),ascending=False).head(28); fig=go.Figure(go.Bar(x=d.values,y=d.index,orientation='h')); return _layout(fig,'Largest ligand rank displacements',520)

def topk_flow_figure(x,scoring,A,B,topk=20):
    _,r,_=_pair_consensus(x,scoring,A,B); a=set(r[A.lower()].nsmallest(topk).index); b=set(r[B.lower()].nsmallest(topk).index); vals=[len(a-b),len(a&b),len(b-a)]
    return _layout(go.Figure(go.Bar(x=['A only','shared','B only'],y=vals)),'Top-k survival',360)

def rrl_curve_figure(x,scoring,A,B,alpha=.05):
    _,_,s=_pair_consensus(x,scoring,A,B); a,b=s[A.lower()].to_numpy(),s[B.lower()].to_numpy(); r=conservative_pair_rrl(a,b,alpha=float(alpha)); vals=[]
    for lab,d in [('A→B',r.A_to_B),('B→A',r.B_to_A)]: vals.append((lab,d.threshold if d.status=='FINITE' else d.support_max))
    return _layout(go.Figure(go.Bar(x=[v[0] for v in vals],y=[v[1] for v in vals])),'Observed-support RRL',360)

def rrl_alpha_sensitivity_figure(x,scoring,A,B):
    _,_,s=_pair_consensus(x,scoring,A,B); a,b=s[A.lower()].to_numpy(),s[B.lower()].to_numpy(); alphas=[.01,.025,.05,.1]; ys=[]
    for z in alphas:
        r=conservative_pair_rrl(a,b,alpha=z); ys.append(r.conservative if r.conservative_status=='FINITE' else r.conservative_lower_bound)
    return _layout(go.Figure(go.Scatter(x=alphas,y=ys,mode='lines+markers')),'RRL criterion sensitivity',360)

def pair_payload(x,scoring,A,B,alpha=.05,topk=20):
    m=pair_metrics(x,scoring,A,B,alpha,topk); comp=[]
    for sc in sorted(x.scoring_function.unique()):
        try:
            z=pair_metrics(x,sc,A,B,alpha,topk); comp.append({'scoring_function':sc,'receptor_rho':z['receptor_spearman_rho'],'seed_rho':z['mean_pair_seed_rho'],'inversion':z['pair_reversal_fraction'],'median_delta_rank':z['median_abs_delta_rank'],'rrl':z['rrl_display']['value'],'rrl_status':z['rrl_display']['status']})
        except Exception: pass
    return {'metrics':m,'rank_figure':_fig_json(rank_transfer_figure(x,scoring,A,B)),'delta_figure':_fig_json(delta_rank_figure(x,scoring,A,B)),'rrl_figure':_fig_json(rrl_curve_figure(x,scoring,A,B,alpha)),'alpha_figure':_fig_json(rrl_alpha_sensitivity_figure(x,scoring,A,B)),'topk_figure':_fig_json(topk_flow_figure(x,scoring,A,B,topk)),'cross_scoring':comp}

def seed_convergence(x,scoring,A,B):
    q=x[(x.scoring_function==scoring.lower()) & x.receptor.isin([A.lower(),B.lower()])]; seeds=sorted(q.seed.unique()); rows=[]
    for n in range(2,len(seeds)+1):
        sub=q[q.seed.isin(seeds[:n])]; rows.append({'n_seeds':n,'receptor_rho':float(analyze_receptor_pair(sub,A,B)['receptor_spearman_rho']),'mean_seed_rho':float(np.mean([seed_reproducibility(sub,r.lower())['mean_seed_rho'] for r in [A,B]]))})
    fig=go.Figure(); fig.add_trace(go.Scatter(x=[r['n_seeds'] for r in rows],y=[r['receptor_rho'] for r in rows],mode='lines+markers',name='receptor')); fig.add_trace(go.Scatter(x=[r['n_seeds'] for r in rows],y=[r['mean_seed_rho'] for r in rows],mode='lines+markers',name='seed'))
    return {'table':rows,'figure':_fig_json(_layout(fig,'Seed-depth convergence',380))}

def validation_panel(x,scoring,A,B,topk=20):
    sub=x[x.scoring_function==scoring.lower()]; cons=consensus_table(sub); scores=cons.pivot(index='ligand_id',columns='receptor',values='score_median'); labels=sub.groupby('ligand_id')['class'].first().reindex(scores.index); y=(labels=='active').astype(int).to_numpy(); m=pair_validity_instability(y,scores[A.lower()].to_numpy(),scores[B.lower()].to_numpy(),topk=int(topk))
    return {'metrics':_json_ready(m),'figure':_fig_json(_layout(go.Figure(go.Bar(x=['A ROC-AUC','B ROC-AUC'],y=[m['roc_auc_A'],m['roc_auc_B']])),'Experimental-label screening',360))}

def screening_landscape(x,scoring,topk=20):
    sub=x[x.scoring_function==scoring.lower()]; rows=[]
    for r in sorted(sub.receptor.unique()):
        q=consensus_table(sub[sub.receptor==r]); labels=sub[sub.receptor==r].groupby('ligand_id')['class'].first().reindex(q.ligand_id); m=screening_metrics((labels=='active').astype(int),q.score_median,topk=int(topk)); rows.append({'receptor':r,**m})
    return {'table':_json_ready(rows),'figure':_fig_json(_layout(go.Figure(go.Bar(x=[r['receptor'] for r in rows],y=[r['roc_auc'] for r in rows])),'Screening landscape',360))}

def ligand_inspector(x,scoring,ligand_id):
    q=x[(x.scoring_function==scoring.lower()) & (x.ligand_id.astype(str)==str(ligand_id))].copy(); return {'table':q.to_dict('records'),'figure':_fig_json(_layout(go.Figure(go.Scatter(x=q.receptor,y=q.score,mode='markers')),'Ligand trace',360))}

def pair_matrix_figure(x,scoring,metric='receptor_spearman_rho'):
    recs=sorted(x.loc[x.scoring_function==scoring.lower(),'receptor'].unique()); mat=np.full((len(recs),len(recs)),np.nan); rows=[]
    for i,a in enumerate(recs):
        mat[i,i]=1
        for j,b in enumerate(recs[i+1:],i+1):
            m=pair_metrics(x,scoring,a,b); v=float(m.get(metric,np.nan)); mat[i,j]=mat[j,i]=v; rows.append({'receptor_a':a,'receptor_b':b,**m})
    fig=go.Figure(go.Heatmap(z=mat,x=recs,y=recs,colorscale='Viridis')); return _layout(fig,f'Pair matrix · {metric}',420),_json_ready(rows)

def seed_vs_receptor_figure(x,scoring):
    rows=[]; recs=sorted(x.loc[x.scoring_function==scoring.lower(),'receptor'].unique())
    for a,b in combinations(recs,2):
        m=pair_metrics(x,scoring,a,b); rows.append(m)
    return _layout(go.Figure(go.Scatter(x=[r['mean_pair_seed_rho'] for r in rows],y=[r['receptor_spearman_rho'] for r in rows],mode='markers')),'Seed repeatability vs receptor transfer',380)

def robustness_panel(x,scoring,A,B,alpha=.05,heldout_repeats=30,subsample_repeats=20):
    _,_,s=_pair_consensus(x,scoring,A,B); a,b=s[A.lower()].to_numpy(),s[B.lower()].to_numpy(); h=heldout_rrl_calibration(a,b,repeats=int(heldout_repeats),alpha=float(alpha)); sizes=[z for z in [100,200,300] if z<=len(a)]; u=subsample_rrl_sensitivity(a,b,sample_sizes=sizes,repeats=int(subsample_repeats),alpha=float(alpha)); return {'heldout':_json_ready(h),'subsample':_json_ready(u)}

def bootstrap_panel(x,scoring,A,B,alpha=.05,n_boot=100,rng_seed=20260813):
    _,_,s=_pair_consensus(x,scoring,A,B); rows=bootstrap_pair_rrl(s[A.lower()].to_numpy(),s[B.lower()].to_numpy(),alpha=float(alpha),n_boot=int(n_boot),rng_seed=int(rng_seed)); q=pd.DataFrame(rows); summary=summarize_bootstrap(rows); finite=q[q.status=='FINITE'] if 'status' in q else pd.DataFrame(); fig=go.Figure(go.Histogram(x=finite.get('rrl_conservative',[])))
    return {'summary':_json_ready(summary),'figure':_fig_json(_layout(fig,f'Ligand-identity bootstrap · B={int(n_boot)}',380)),'table':q.to_dict('records')}

def demo_atlas(master,biological):
    x=master.copy(); fig=go.Figure()
    for sc,q in x.groupby('scoring_function'): fig.add_trace(go.Scatter(x=q.mean_within_seed_rho,y=q.receptor_spearman_rho,mode='markers',name=str(sc)))
    return {'figure':_fig_json(_layout(fig,'Canonical DockRIFT project atlas',430)),'summary':{'records':int(len(x)),'targets':int(x.biological_target.nunique())},'target_cards':[],'table':x.to_dict('records')}

def provenance_payload(session):
    import numpy,pandas,scipy,sklearn,plotly
    return {'dockrift_version':__version__,'source_name':session.source_name,'source_sha256':session.sha256,'source_kind':session.source_kind,'python':sys.version.split()[0],'platform':platform.platform(),'numpy':numpy.__version__,'pandas':pandas.__version__,'scipy':scipy.__version__,'scikit_learn':sklearn.__version__,'plotly':plotly.__version__,'method':{'ranking':'lower docking score is better; average ranks for ties','consensus':'median score across seeds per ligand/receptor','rrl':'directional observed-support isotonic reversal model; no extrapolation; conservative pair propagation','interpretation':'RRL is a docking-score-scale rank-resolution descriptor, not a binding-free-energy uncertainty'}}

def export_bundle(session,scoring,A,B,alpha=.05,topk=20):
    out=Path(tempfile.mkdtemp(prefix='dockrift_studio_')); zpath=out/f'DockRIFT_{A.upper()}_{B.upper()}_{scoring}_bundle.zip'; p=pair_payload(session.data,scoring,A,B,alpha,topk); files=[]
    n=out/'normalized_scores.csv'; session.data.to_csv(n,index=False); files.append(n)
    m=out/'pair_metrics.csv'; pd.DataFrame([p['metrics']]).to_csv(m,index=False); files.append(m)
    pr=out/'provenance.json'; pr.write_text(json.dumps(provenance_payload(session),indent=2)); files.append(pr)
    mc=out/'methods_capsule.txt'; mc.write_text('DockRIFT Studio v'+__version__+'\nRRL: directional observed-support isotonic reversal model; no extrapolation.\n'); files.append(mc)
    h=out/'interactive_report.html'; h.write_text('<!doctype html><meta charset="utf-8"><title>DockRIFT report</title><h1>DockRIFT Studio report</h1>'); files.append(h)
    with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED) as z:
        for f in files: z.write(f,f.name)
    return str(zpath)
