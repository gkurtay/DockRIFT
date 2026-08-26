from __future__ import annotations

import argparse
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
import json
import mimetypes
from pathlib import Path
import secrets
import threading
import urllib.parse
import webbrowser

import plotly.offline as po

from .logic import (
    SessionData,
    bootstrap_panel,
    data_health,
    dataset_profile,
    demo_atlas,
    export_bundle,
    ligand_inspector,
    load_demo,
    pair_matrix_figure,
    pair_payload,
    provenance_payload,
    read_score_bytes,
    robustness_panel,
    screening_landscape,
    seed_convergence,
    seed_vs_receptor_figure,
    validation_panel,
    _fig_json,
)
from .. import __version__

STATIC = Path(__file__).resolve().parent / "static"
_LOCK = threading.RLock()
_SESSION, _DEMO_MASTER, _DEMO_BIO = load_demo()
_DOWNLOADS: dict[str, Path] = {}


def _safe_json(obj):
    return json.dumps(obj, allow_nan=False).encode("utf-8")


def _session_snapshot():
    with _LOCK:
        return _SESSION


def _set_session(s: SessionData):
    global _SESSION
    with _LOCK:
        _SESSION = s


class Handler(BaseHTTPRequestHandler):
    server_version = f"DockRIFTStudio/{__version__}"

    def log_message(self, fmt, *args):
        print("[DockRIFT Studio] " + fmt % args)

    def _send(self, code: int, content: bytes, content_type="application/json; charset=utf-8", headers=None):
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(content)))
        self.send_header("Cache-Control", "no-store" if content_type.startswith("application/json") else "public, max-age=3600")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("X-Frame-Options", "SAMEORIGIN")
        self.send_header("Referrer-Policy", "no-referrer")
        if headers:
            for k,v in headers.items(): self.send_header(k,v)
        self.end_headers()
        self.wfile.write(content)

    def _json(self, payload, code=200):
        try:
            self._send(code, _safe_json(payload))
        except Exception as e:
            self._send(500, _safe_json({"ok":False,"error":f"serialization failure: {e}"}))

    def _body_json(self):
        length = int(self.headers.get("Content-Length", "0"))
        if length > 100 * 1024 * 1024:
            raise ValueError("Request exceeds the 100 MB local Studio limit")
        body = self.rfile.read(length) if length else b"{}"
        return json.loads(body.decode("utf-8"))

    def do_GET(self):
        u = urllib.parse.urlparse(self.path)
        path = u.path
        if path == "/api/status":
            s=_session_snapshot()
            return self._json({"ok":True,"version":__version__,"profile":dataset_profile(s.data,s.source_name,s.sha256),"local_only":True})
        if path == "/api/demo-atlas":
            return self._json({"ok":True,**demo_atlas(_DEMO_MASTER,_DEMO_BIO)})
        if path == "/assets/plotly.min.js":
            return self._send(200, po.get_plotlyjs().encode("utf-8"), "application/javascript; charset=utf-8")
        if path.startswith("/download/"):
            token=path.rsplit("/",1)[-1]
            p=_DOWNLOADS.get(token)
            if not p or not p.is_file(): return self._json({"ok":False,"error":"Download token expired or invalid"},404)
            data=p.read_bytes()
            return self._send(200,data,"application/zip",{"Content-Disposition":f'attachment; filename="{p.name}"'})
        if path == "/": path="/index.html"
        rel=Path(path.lstrip("/"))
        target=(STATIC/rel).resolve()
        try: target.relative_to(STATIC.resolve())
        except ValueError: return self._json({"ok":False,"error":"Forbidden"},403)
        if not target.is_file(): return self._json({"ok":False,"error":"Not found"},404)
        ctype=mimetypes.guess_type(str(target))[0] or "application/octet-stream"
        return self._send(200,target.read_bytes(),ctype)

    def do_POST(self):
        try:
            payload=self._body_json()
            path=urllib.parse.urlparse(self.path).path
            if path == "/api/demo":
                s,_,_=load_demo();_set_session(s)
                return self._json({"ok":True,"profile":dataset_profile(s.data,s.source_name,s.sha256),"health":data_health(s.data)})
            if path == "/api/upload":
                filename=str(payload.get("filename") or "scores.csv")
                text=payload.get("content")
                if text is None: raise ValueError("No file content received")
                s=read_score_bytes(filename,text.encode("utf-8"));_set_session(s)
                return self._json({"ok":True,"profile":dataset_profile(s.data,s.source_name,s.sha256),"health":data_health(s.data)})
            s=_session_snapshot();x=s.data
            if path == "/api/profile":
                return self._json({"ok":True,"profile":dataset_profile(x,s.source_name,s.sha256),"health":data_health(x)})
            if path == "/api/health":
                return self._json({"ok":True,**data_health(x)})
            if path == "/api/pair":
                sc=str(payload["scoring"]);A=str(payload["receptor_a"]);B=str(payload["receptor_b"]);alpha=float(payload.get("alpha",.05));topk=int(payload.get("topk",20))
                return self._json({"ok":True,**pair_payload(x,sc,A,B,alpha,topk)})
            if path == "/api/atlas":
                sc=str(payload["scoring"]);metric=str(payload.get("metric","receptor_spearman_rho"));fig,table=pair_matrix_figure(x,sc,metric)
                return self._json({"ok":True,"matrix_figure":_fig_json(fig),"transfer_figure":_fig_json(seed_vs_receptor_figure(x,sc)),"table":table})
            if path == "/api/validation":
                return self._json({"ok":True,**validation_panel(x,str(payload["scoring"]),str(payload["receptor_a"]),str(payload["receptor_b"]),int(payload.get("topk",20)))})
            if path == "/api/screening-landscape":
                return self._json({"ok":True,**screening_landscape(x,str(payload["scoring"]),int(payload.get("topk",20)))})
            if path == "/api/seed-convergence":
                return self._json({"ok":True,**seed_convergence(x,str(payload["scoring"]),str(payload["receptor_a"]),str(payload["receptor_b"]))})
            if path == "/api/ligand":
                return self._json({"ok":True,**ligand_inspector(x,str(payload["scoring"]),str(payload["ligand_id"]))})
            if path == "/api/robustness":
                return self._json({"ok":True,**robustness_panel(x,str(payload["scoring"]),str(payload["receptor_a"]),str(payload["receptor_b"]),float(payload.get("alpha",.05)),int(payload.get("heldout_repeats",30)),int(payload.get("subsample_repeats",20)))})
            if path == "/api/bootstrap":
                return self._json({"ok":True,**bootstrap_panel(x,str(payload["scoring"]),str(payload["receptor_a"]),str(payload["receptor_b"]),float(payload.get("alpha",.05)),int(payload.get("n_boot",100)),int(payload.get("rng_seed",20260813)))})
            if path == "/api/provenance":
                return self._json({"ok":True,"provenance":provenance_payload(s)})
            if path == "/api/export":
                p=Path(export_bundle(s,str(payload["scoring"]),str(payload["receptor_a"]),str(payload["receptor_b"]),float(payload.get("alpha",.05)),int(payload.get("topk",20))))
                token=secrets.token_urlsafe(16);_DOWNLOADS[token]=p
                return self._json({"ok":True,"download_url":f"/download/{token}","filename":p.name})
            return self._json({"ok":False,"error":"Unknown API endpoint"},404)
        except Exception as e:
            return self._json({"ok":False,"error":f"{type(e).__name__}: {e}"},400)


def serve(host="127.0.0.1", port=8765, open_browser=True):
    server=ThreadingHTTPServer((host,int(port)),Handler)
    url=f"http://{host}:{int(port)}/"
    print(f"DockRIFT Studio {__version__}")
    print(f"Local interface: {url}")
    print("Local-first mode: uploaded files remain in this Python process.")
    if open_browser:
        threading.Timer(.7,lambda:webbrowser.open(url)).start()
    try: server.serve_forever()
    except KeyboardInterrupt: pass
    finally: server.server_close()


def main(argv=None):
    p=argparse.ArgumentParser(description="Launch the local DockRIFT Studio browser interface")
    p.add_argument("--host",default="127.0.0.1")
    p.add_argument("--port",type=int,default=8765)
    p.add_argument("--no-browser",action="store_true")
    args=p.parse_args(argv)
    serve(args.host,args.port,not args.no_browser)


if __name__ == "__main__": main()
