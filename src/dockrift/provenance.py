from __future__ import annotations
from pathlib import Path
import hashlib


def sha256_file(path:str|Path)->str:
    p=Path(path);h=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()


def write_sha256_manifest(paths,out_path:str|Path)->Path:
    out=Path(out_path);out.parent.mkdir(parents=True,exist_ok=True);lines=[]
    for item in paths:
        p=Path(item).resolve()
        if not p.is_file():raise FileNotFoundError(p)
        lines.append(f"{sha256_file(p)}  {p}")
    out.write_text("\n".join(lines)+"\n",encoding="utf-8");return out


def verify_sha256_manifest(path:str|Path)->list[dict]:
    p=Path(path);rows=[]
    for line in p.read_text(errors="replace").splitlines():
        if not line.strip():continue
        digest,name=line.split(None,1);name=name.strip();f=Path(name);actual=sha256_file(f) if f.is_file() else None;rows.append({"path":name,"expected":digest,"actual":actual,"ok":actual==digest})
    return rows
