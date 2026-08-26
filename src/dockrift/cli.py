from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import sys

from . import __version__
from .report import write_analysis


def build_parser():
    p=argparse.ArgumentParser(prog="dockrift",description="DockRIFT rank-stability / transferability toolkit")
    p.add_argument("--version",action="version",version=f"DockRIFT {__version__}")
    sp=p.add_subparsers(dest="command")
    a=sp.add_parser("analyze",help="Analyze a balanced long-form score table")
    a.add_argument("scores");a.add_argument("--out",required=True);a.add_argument("--alpha",type=float,default=.05);a.add_argument("--scoring-function",default=None)
    s=sp.add_parser("studio",help="Launch the custom local DockRIFT Studio browser GUI")
    s.add_argument("--host",default="127.0.0.1");s.add_argument("--port",type=int,default=8765);s.add_argument("--no-browser",action="store_true")
    e=sp.add_parser("extract",help="Extract mode-1 scores from a DockRIFT/Vina task manifest")
    e.add_argument("manifest");e.add_argument("--out",required=True)
    f=sp.add_parser("freeze",help="Write a SHA256 provenance manifest for exact files")
    f.add_argument("files",nargs="+");f.add_argument("--out",required=True)
    v=sp.add_parser("verify",help="Verify a DockRIFT SHA256 provenance manifest")
    v.add_argument("manifest")
    sp.add_parser("doctor",help="Print environment and GUI readiness checks")
    return p


def doctor():
    import numpy,pandas,scipy,sklearn,plotly
    from .gui.server import STATIC
    checks={
        "dockrift":__version__,"python":sys.version.split()[0],"platform":platform.platform(),
        "numpy":numpy.__version__,"pandas":pandas.__version__,"scipy":scipy.__version__,"scikit_learn":sklearn.__version__,"plotly":plotly.__version__,
        "studio_static":all((STATIC/f).is_file() for f in ["index.html","style.css","app.js"]),
        "demo_bundle": False,
    }
    try:
        from .gui.logic import load_demo
        demo, master, bio = load_demo()
        checks["demo_bundle"] = bool(len(demo.data) == 4800 and len(master) == 60 and len(bio) == 14)
    except Exception:
        checks["demo_bundle"] = False
    checks["ready"]=bool(checks["studio_static"] and checks["demo_bundle"])
    print(json.dumps(checks,indent=2))
    return checks


def main(argv=None):
    args=build_parser().parse_args(argv)
    if not args.command:
        build_parser().print_help();return 0
    if args.command=="analyze":
        df=write_analysis(args.scores,args.out,alpha=args.alpha,scoring_function=args.scoring_function);print(f"DockRIFT analysis complete: {len(df)} receptor pairs -> {args.out}")
    elif args.command=="studio":
        from .gui.server import serve
        serve(args.host,args.port,not args.no_browser)
    elif args.command=="extract":
        from .extract import score_table_from_manifest
        df=score_table_from_manifest(args.manifest,args.out);print(f"Extracted {len(df)} mode-1 scores -> {args.out}")
    elif args.command=="freeze":
        from .provenance import write_sha256_manifest
        p=write_sha256_manifest(args.files,args.out);print(f"Wrote SHA256 manifest: {p}")
    elif args.command=="verify":
        from .provenance import verify_sha256_manifest
        rows=verify_sha256_manifest(args.manifest);
        for r in rows:print(f"{r['path']}: {'OK' if r['ok'] else 'FAIL'}")
        if not all(r['ok'] for r in rows):return 2
    elif args.command=="doctor": doctor()
    return 0

if __name__=="__main__": raise SystemExit(main())
