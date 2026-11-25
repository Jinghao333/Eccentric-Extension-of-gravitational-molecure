"""Local read-only browser dashboard for Xu Appendix-A v4 runs.

Usage (from fig8_minimal/):
  python -m analysis_dashboard.server --runs runs --port 8765
Open http://127.0.0.1:8765

Does NOT write scientific products. Scientific API is in core.py and is
available unchanged to notebooks, tests, and independent analysis scripts.
"""
from __future__ import annotations
import argparse
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
import json
from pathlib import Path
from urllib.parse import urlparse,parse_qs
import traceback
import numpy as np
from .core import load_run, run_overview,channel_atlas,density_slice,channel_detail

STATIC=Path(__file__).parent/'static'


def run_folders(root):
    root=Path(root)
    if (root/'config.json').is_file() and (root/'physics.h5').is_file():
        return {root.name:root.resolve()}
    return {p.name:p.resolve() for p in sorted(root.glob('*')) if p.is_dir() and
            (p/'config.json').is_file() and (p/'physics.h5').is_file()}


def create_handler(root):
    root=Path(root).resolve()
    class Handler(BaseHTTPRequestHandler):
        def send_bytes(self,body,kind='application/json',status=200):
            self.send_response(status)
            self.send_header('Content-Type',kind)
            self.send_header('Content-Length',str(len(body)))
            self.send_header('Cache-Control','no-store')
            self.send_header('X-Content-Type-Options','nosniff')
            self.end_headers();self.wfile.write(body)

        def send_json(self,obj,status=200):
            def convert(x):
                if isinstance(x,(np.integer,np.floating)):return x.item()
                if isinstance(x,np.ndarray):return x.tolist()
                raise TypeError(type(x).__name__)
            body=json.dumps(obj,default=convert,allow_nan=False,separators=(',',':')).encode('utf-8')
            self.send_bytes(body,status=status)

        def do_GET(self):
            url=urlparse(self.path);path=url.path;qs=parse_qs(url.query)
            try:
                if path=='/' or path=='/index.html':
                    self.send_bytes((STATIC/'index.html').read_bytes(),'text/html; charset=utf-8');return
                if path=='/style.css':
                    self.send_bytes((STATIC/'style.css').read_bytes(),'text/css; charset=utf-8');return
                if path=='/app.js':
                    self.send_bytes((STATIC/'app.js').read_bytes(),'text/javascript; charset=utf-8');return
                if path=='/plotly.min.js':
                    from plotly.offline import get_plotlyjs
                    self.send_bytes(get_plotlyjs().encode('utf-8'),'text/javascript; charset=utf-8');return
                if path=='/api/runs':
                    self.send_json({'runs':list(run_folders(root))});return
                folders=run_folders(root)
                name=qs.get('run',[''])[0]
                if name not in folders:
                    self.send_json({'error':'Unknown or incomplete run. Choose a completed run in the selector.'},404);return
                result=load_run(folders[name])
                if path=='/api/overview':
                    self.send_json(run_overview(result));return
                if path=='/api/atlas':
                    self.send_json(channel_atlas(result));return
                if path=='/api/geometry':
                    n=max(90,min(300,int(qs.get('n',['170'])[0])))
                    ext=qs.get('extent',[None])[0]
                    self.send_json(density_slice(result['run']['config'],n=n,
                                                 extent=float(ext) if ext is not None else None));return
                if path=='/api/channel':
                    ell=int(qs['ell'][0]);mf=int(qs['m'][0]);g=int(qs['g'][0])
                    self.send_json(channel_detail(result,ell,mf,g,True));return
                self.send_json({'error':'Unknown endpoint'},404)
            except (KeyError,ValueError,ArithmeticError) as exc:
                self.send_json({'error':str(exc)},422)
            except BrokenPipeError:
                pass
            except Exception as exc:
                traceback.print_exc()
                self.send_json({'error':f'{type(exc).__name__}: {exc}'},500)

        def log_message(self,format,*args):
            if not self.path.startswith('/api/'):
                return
            super().log_message(format,*args)
    return Handler


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs',type=Path,default=Path('runs'))
    parser.add_argument('--host',default='127.0.0.1')
    parser.add_argument('--port',type=int,default=8765)
    args=parser.parse_args()
    if not args.runs.exists():
        raise SystemExit(f'Run folder {args.runs!s} does not exist')
    server=ThreadingHTTPServer((args.host,args.port),create_handler(args.runs))
    print(f'Xu Appendix-A dashboard at http://{args.host}:{args.port}/')
    try:server.serve_forever()
    except KeyboardInterrupt:pass
    finally:server.server_close()

if __name__=='__main__':main()
