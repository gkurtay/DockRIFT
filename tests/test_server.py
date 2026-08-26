import json,threading,urllib.request
from http.server import ThreadingHTTPServer
from dockrift.gui.server import Handler

def test_local_server_status_and_static():
    srv=ThreadingHTTPServer(("127.0.0.1",0),Handler);port=srv.server_address[1]
    th=threading.Thread(target=srv.serve_forever,daemon=True);th.start()
    try:
        with urllib.request.urlopen(f"http://127.0.0.1:{port}/api/status") as r:
            j=json.loads(r.read());assert j["ok"] and j["profile"]["ligands"]==400
        with urllib.request.urlopen(f"http://127.0.0.1:{port}/") as r:
            html=r.read().decode();assert "DockRIFT Studio" in html and "Pair Explorer" in html
        with urllib.request.urlopen(f"http://127.0.0.1:{port}/assets/plotly.min.js") as r:
            assert len(r.read())>1000000
    finally:
        srv.shutdown();srv.server_close()
