from __future__ import annotations
import argparse, json, mimetypes, secrets, threading, urllib.parse, webbrowser
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from pathlib import Path
import plotly.offline as po
from .logic import (load_demo,dataset_profile,data_health,read_score_bytes,pair_payload,
    validation_panel,seed_convergence,provenance_payload,export_bundle,bootstrap_panel,
    robustness_panel,pair_matrix_figure,seed_vs_receptor_figure,screening_landscape,
    ligand_inspector,demo_atlas,_fig_json)
from .. import __version__
STATIC=Path(__file__).resolve().parent/'static';_LOCK=threading.RLock();_SESSION,_DEMO_MASTER,_DEMO_BIO=load_demo();_DOWNLOADS={}
def _session():
    with _LOCK:return _SESSION
def _set(s):
    global _SESSION
    with _LOCK:_SESSION=s
class Handler(BaseHTTPRequestHandler):
    server_version=f'DockRIFTStudio/{__version__}'
    def log_message(self,fmt,*args):print('[DockRIFT Studio] '+fmt%args)
    def _send(self,code,data,ctype='application/json; charset=utf-8',headers=None):
        self.send_response(code);self.send_header('Content-Type',ctype);self.send_header('Content-Length',str(len(data)));self.send_header('Cache-Control','no-store' if ctype.startswith('application/json') else 'public, max-age=3600');self.send_header('X-Content-Type-Options','nosniff');self.send_header('X-Frame-Options','SAMEORIGIN');self.send_header('Referrer-Policy','no-referrer')
        if headers:
            for k,v in headers.items():self.send_header(k,v)
        self.end_headers();self.wfile.write(data)
    def _json(self,obj,code=200):self._send(code,json.dumps(obj,allow_nan=False).encode())
    def _body(self):
        n=int(self.headers.get('Content-Length','0'))
        if n>100*1024*1024:raise ValueError('Request exceeds the 100 MB local Studio limit')
        return json.loads(self.rfile.read(n).decode() if n else '{}')
    def do_GET(self):
        path=urllib.parse.urlparse(self.path).path
        if path=='/api/status':
            s=_session();return self._json({'ok':True,'version':__version__,'profile':dataset_profile(s.data,s.source_name,s.sha256),'local_only':True})
        if path=='/api/demo-atlas':return self._json({'ok':True,**demo_atlas(_DEMO_MASTER,_DEMO_BIO)})
        if path=='/assets/plotly.min.js':return self._send(200,po.get_plotlyjs().encode(),'application/javascript; charset=utf-8')
        if path.startswith('/download/'):
            p=_DOWNLOADS.get(path.rsplit('/',1)[-1])
            if not p or not Path(p).is_file():return self._json({'ok':False,'error':'invalid or expired download'},404)
            b=Path(p).read_bytes();return self._send(200,b,'application/zip',{'Content-Disposition':f'attachment; filename="{Path(p).name}"'})
        if path=='/':path='/index.html'
        target=(STATIC/path.lstrip('/')).resolve()
        try:target.relative_to(STATIC.resolve())
        except ValueError:return self._json({'ok':False,'error':'forbidden'},403)
        if not target.is_file():return self._json({'ok':False,'error':'not found'},404)
        return self._send(200,target.read_bytes(),mimetypes.guess_type(str(target))[0] or 'application/octet-stream')
    def do_POST(self):
        try:
            p=self._body();path=urllib.parse.urlparse(self.path).path
            if path=='/api/demo':
                s,_,_=load_demo();_set(s);return self._json({'ok':True,'profile':dataset_profile(s.data,s.source_name,s.sha256),'health':data_health(s.data)})
            if path=='/api/upload':
                s=read_score_bytes(str(p.get('filename','scores.csv')),str(p['content']).encode());_set(s);return self._json({'ok':True,'profile':dataset_profile(s.data,s.source_name,s.sha256),'health':data_health(s.data)})
            s=_session();x=s.data
            if path=='/api/profile':return self._json({'ok':True,'profile':dataset_profile(x,s.source_name,s.sha256),'health':data_health(x)})
            if path=='/api/health':return self._json({'ok':True,**data_health(x)})
            if path=='/api/pair':return self._json({'ok':True,**pair_payload(x,p['scoring'],p['receptor_a'],p['receptor_b'],float(p.get('alpha',.05)),int(p.get('topk',20)))})
            if path=='/api/atlas':
                fig,table=pair_matrix_figure(x,p['scoring'],p.get('metric','receptor_spearman_rho'));return self._json({'ok':True,'matrix_figure':_fig_json(fig),'transfer_figure':_fig_json(seed_vs_receptor_figure(x,p['scoring'])),'table':table})
            if path=='/api/validation':return self._json({'ok':True,**validation_panel(x,p['scoring'],p['receptor_a'],p['receptor_b'],int(p.get('topk',20)))})
            if path=='/api/screening-landscape':return self._json({'ok':True,**screening_landscape(x,p['scoring'],int(p.get('topk',20)))})
            if path=='/api/seed-convergence':return self._json({'ok':True,**seed_convergence(x,p['scoring'],p['receptor_a'],p['receptor_b'])})
            if path=='/api/ligand':return self._json({'ok':True,**ligand_inspector(x,p['scoring'],p['ligand_id'])})
            if path=='/api/robustness':return self._json({'ok':True,**robustness_panel(x,p['scoring'],p['receptor_a'],p['receptor_b'],float(p.get('alpha',.05)),int(p.get('heldout_repeats',30)),int(p.get('subsample_repeats',20)))})
            if path=='/api/bootstrap':return self._json({'ok':True,**bootstrap_panel(x,p['scoring'],p['receptor_a'],p['receptor_b'],float(p.get('alpha',.05)),int(p.get('n_boot',20)),int(p.get('rng_seed',20260813)))})
            if path=='/api/provenance':return self._json({'ok':True,**provenance_payload(s)})
            if path=='/api/export':
                token=secrets.token_urlsafe(12);_DOWNLOADS[token]=Path(export_bundle(s,p['scoring'],p['receptor_a'],p['receptor_b'],float(p.get('alpha',.05)),int(p.get('topk',20))));return self._json({'ok':True,'download_url':'/download/'+token})
            return self._json({'ok':False,'error':'unknown endpoint'},404)
        except Exception as e:return self._json({'ok':False,'error':str(e)},400)
def serve(host='127.0.0.1',port=8765,open_browser=True):
    srv=ThreadingHTTPServer((host,int(port)),Handler);url=f'http://{host}:{port}/';print(f'DockRIFT Studio {__version__}: {url}')
    if open_browser:threading.Timer(.35,lambda:webbrowser.open(url)).start()
    try:srv.serve_forever()
    except KeyboardInterrupt:pass
    finally:srv.server_close()
def main(argv=None):
    ap=argparse.ArgumentParser();ap.add_argument('--host',default='127.0.0.1');ap.add_argument('--port',type=int,default=8765);ap.add_argument('--no-browser',action='store_true');a=ap.parse_args(argv);serve(a.host,a.port,not a.no_browser)
if __name__=='__main__':main()
