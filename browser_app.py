#!/usr/bin/env python3
"""Private, local browser interface for PlasmidTrace (standard library only)."""
from __future__ import annotations
import argparse, json, os, secrets, shlex, shutil, subprocess, sys, threading, time, webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit
from local_app import Backend, CONDA_INIT, Job, ROOT, make_archive, required_tools
from automatic_setup import SetupTask, apply_runtime, validate_screening

ASSETS=ROOT/'app_assets'
BOOLS={'run_barrnap','run_mob','run_plasmidfinder','automatic','exploratory'}
DEFAULTS=dict(automatic=True,exploratory=False,runtime_python='',workflow='screen',assembly='',reads1='',reads2='',previous='',output_root=str(ROOT/'results'),threads='2',memory='7',backend='local',distro='Ubuntu-24.04',environment='',recovery_python='python',mob_executable='mob_recon',plasmidfinder_executable='plasmidfinder.py',plasmidfinder_python='',run_barrnap=False,run_mob=False,run_plasmidfinder=False)
for name in ('rrna_gff','mob_report','gfa','annotation_genbank','annotation_gff','annotation_fasta','plasmidfinder_results','plasmidfinder_db','mob_db'):DEFAULTS[name]=''
for tool in ('amrfinder','virulencefinder','geneminer'):
    DEFAULTS[tool+'_mode']='off';DEFAULTS[tool+'_result']=''
DEFAULTS.update(amrfinder_exe='amrfinder',amrfinder_db='',virulencefinder_python='python',virulencefinder_db='',virulencefinder_blast='blastn',geneminer_genbank='',geneminer_workbook='',geneminer_python='python',geneminer_script='',geneminer_db='',geneminer_metadata='',geneminer_rules='')

def clean_settings(data):
    if not isinstance(data,dict):raise ValueError('Settings must be an object')
    result=dict(DEFAULTS)
    for k,v in data.items():
        if k not in DEFAULTS:raise ValueError('Unknown setting: '+k)
        if k in BOOLS:
            if not isinstance(v,bool):raise ValueError(k+' must be true or false')
        elif not isinstance(v,str) or len(v)>8192:raise ValueError('Invalid value for '+k)
        result[k]=v.strip() if isinstance(v,str) else v
    if result['backend'] not in ('local','wsl'):raise ValueError('Choose a run environment')
    if result['workflow'] not in ('screen','reconstruct','reuse','recover','cargo'):raise ValueError('Choose a workflow')
    if not result['output_root']:raise ValueError('Choose a result folder')
    return result

def run_settings(saved):
    s=dict(saved)
    if s['workflow']=='screen':
        for tool in ('amrfinder','virulencefinder','geneminer'):s[tool+'_mode']='off'
    if s['workflow']=='cargo':
        for k in ('run_barrnap','run_mob','run_plasmidfinder'):s[k]=False
        for k in ('rrna_gff','mob_report','plasmidfinder_results','annotation_genbank','annotation_gff','annotation_fasta'):s[k]=''
    return s

def read_tail(path,limit=60000):
    try:
        with path.open('rb') as f:
            f.seek(max(0,path.stat().st_size-limit));return f.read().decode('utf-8',errors='replace')
    except OSError:return ''

class State:
    def __init__(self,root=ROOT):
        self.root=Path(root);self.settings_file=self.root/'app_settings.json';self.jobs={};self.lock=threading.RLock();self.picker_lock=threading.Lock();self.zip_lock=threading.Lock();self.active=None;self.setup_task=None;self.runtime_result={};self.setup_distro=None
        try:self.saved=clean_settings(json.loads(self.settings_file.read_text(encoding='utf-8')))
        except (OSError,ValueError):self.saved=dict(DEFAULTS,output_root=str(self.root/'results'))
    def start(self,data):
        saved=clean_settings(data);s=run_settings(saved)
        if s['automatic'] and s['workflow']!='cargo':
            if s['distro']!=self.setup_distro:raise ValueError('Run automatic detection for the selected Linux distribution first.')
            s=apply_runtime(s,self.runtime_result)
        validate_screening(s)
        with self.lock:
            if self.setup_task and not self.setup_task.finished.is_set():raise ValueError('Let setup finish before starting an analysis.')
            if self.active and not self.jobs[self.active].finished.is_set():raise ValueError('A run is already active. Wait or cancel it first.')
            backend=Backend(s['backend'],s['distro'],s['environment'],python=s.get('runtime_python') or None)
            job=Job(s,backend,s['output_root']);key=job.folder.name
            self.settings_file.write_text(json.dumps(saved,indent=2),encoding='utf-8');self.saved=saved
            self.jobs[key]=job;self.active=key
            threading.Thread(target=job.run,daemon=True).start()
            return key
    def status(self,key):
        job=self.jobs[key];done=job.finished.is_set()
        logs=read_tail(job.folder/'console.log')
        try:children=sorted((job.folder/'output').rglob('*.log'),key=lambda p:p.stat().st_mtime)[-3:]
        except OSError:children=[]
        for p in children:logs+='\n\n— '+str(p.relative_to(job.folder))+' —\n'+read_tail(p,6000)
        assessment={}
        try:assessment=json.loads((job.folder/'assessment.json').read_text())
        except (OSError,ValueError):pass
        return dict(assessment=assessment,id=key,state=job.state,done=done,log=logs,folder=str(job.folder),workflow=job.settings['workflow'],report=done and (job.folder/'report.html').is_file(),excel=done and (job.folder/'PlasmidTrace_results.xlsx').is_file(),elapsed=max(0,int(time.time()-job.folder.stat().st_mtime)))
    def setup(self,data,install=False):
        s=clean_settings(data)
        with self.lock:
            if self.active and not self.jobs[self.active].finished.is_set():raise ValueError('An analysis is active. Wait or cancel it before changing setup.')
            if self.setup_task and not self.setup_task.finished.is_set():raise ValueError('Setup is already running.')
            self.setup_task=SetupTask(s['distro'],self.root,install=install)
            # A changed distribution invalidates any earlier detection result.
            self.runtime_result={}
            self.setup_distro=s['distro']
            threading.Thread(target=self.setup_task.run,daemon=True).start()
            return self.setup_task.id
    def setup_status(self):
        if not self.setup_task:return {'state':'not_checked','done':True,'result':{}}
        d=self.setup_task.status()
        if d['done'] and d['state']=='finished' and d['result'].get('ready'):self.runtime_result=d['result']
        return d
    def check(self,data):
        s=run_settings(clean_settings(data))
        if s['automatic'] and s['workflow']!='cargo':
            if s['distro']!=self.setup_distro:raise ValueError('Run automatic detection for the selected Linux distribution first.')
            s=apply_runtime(s,self.runtime_result)
        b=Backend(s['backend'],s['distro'],s['environment'],python=s.get('runtime_python') or None);info=b.probe();required=required_tools(s)
        extra=[b.path(n) for n in required if b.path(n) not in info['tools']]
        if extra:info['tools'].update(b.probe(extra)['tools'])
        missing=[n for n in required if not info['tools'].get(b.path(n))]
        return dict(**info,missing=missing,ready=not missing,backend=b.kind)
    def discover(self,data):
        s=clean_settings(data);b=Backend('wsl',s['distro'])
        r=subprocess.run(b.shell(CONDA_INIT+'conda env list --json'),capture_output=True,text=True,timeout=45)
        if r.returncode:raise ValueError(r.stderr.strip() or 'Conda was not found in this WSL distribution.')
        envs=json.loads(r.stdout)['envs']
        # Inspect installed executables without activating every environment.
        names=required_tools(run_settings(s)) or ['barrnap','mob_recon','unicycler','blastn','bwa','samtools']
        code='import os,json; envs='+repr(envs)+'; names='+repr([Path(n).name for n in names])+'; scores={e:sum(os.access(os.path.join(e,"bin",n),os.X_OK) for n in names) for e in envs}; print(json.dumps(scores))'
        scores={}
        try:
            check=subprocess.run(b.shell('python3 -c '+shlex.quote(code)),capture_output=True,text=True,timeout=15)
            if check.returncode==0:scores=json.loads(check.stdout)
        except (ValueError,OSError,subprocess.SubprocessError):pass
        recommended=max(scores,key=scores.get) if scores and max(scores.values())>0 else ''
        return {'environments':envs,'recommended':recommended}
    def pick(self,kind):
        if kind not in ('file','folder'):raise ValueError('Invalid selection type')
        if not self.picker_lock.acquire(blocking=False):raise ValueError('A file-selection window is already open. Check your taskbar.')
        try:
            kwargs={'creationflags':subprocess.CREATE_NO_WINDOW} if os.name=='nt' else {}
            r=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--pick',kind],capture_output=True,text=True,timeout=600,**kwargs)
            if r.returncode:raise ValueError('The file picker could not open. You can paste the full path into the field instead.')
            return {'path':json.loads(r.stdout)}
        finally:self.picker_lock.release()

class Handler(BaseHTTPRequestHandler):
    def log_message(self,*args):pass
    def route(self):
        if self.headers.get('Host')!=f'127.0.0.1:{self.server.server_port}':raise PermissionError('Invalid local host')
        path=urlsplit(self.path).path;prefix='/'+self.server.token+'/'
        if not path.startswith(prefix):raise PermissionError('Use the private link printed in the launcher')
        return path[len(prefix):]
    def headers_for(self,length,kind,status=200,attachment=None,report=False):
        self.send_response(status);self.send_header('Content-Type',kind);self.send_header('Content-Length',str(length));self.send_header('Cache-Control','no-store');self.send_header('X-Content-Type-Options','nosniff');self.send_header('Referrer-Policy','no-referrer')
        self.send_header('Content-Security-Policy',"sandbox; default-src 'none'; style-src 'unsafe-inline'; frame-ancestors 'none'" if report else "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; frame-ancestors 'none'; base-uri 'none'")
        if attachment:self.send_header('Content-Disposition','attachment; filename="'+attachment+'"')
        self.end_headers()
    def send(self,data,status=200,kind='application/json'):
        if kind=='application/json':data=json.dumps(data).encode()
        elif isinstance(data,str):data=data.encode()
        self.headers_for(len(data),kind,status);self.wfile.write(data)
    def file(self,p,kind,attachment=None,report=False):
        with p.open('rb') as f:
            self.headers_for(p.stat().st_size,kind,attachment=attachment,report=report);shutil.copyfileobj(f,self.wfile)
    def do_GET(self):
        try:
            path=self.route();state=self.server.state
            if path=='':return self.file(ASSETS/'index.html','text/html; charset=utf-8')
            if path in ('app.css','app.js'):return self.file(ASSETS/path,'text/css; charset=utf-8' if path.endswith('css') else 'text/javascript; charset=utf-8')
            if path=='api/settings':return self.send({'settings':state.saved,'active':state.active,'version':'0.3.0-preview'})
            if path=='api/setup/status':return self.send(state.setup_status())
            if path.startswith('api/status/'):return self.send(state.status(path.split('/')[-1]))
            if path.startswith(('report/','excel/','download/')):
                job=state.jobs[path.split('/')[-1]]
                if not job.finished.is_set():raise ValueError('Wait until the run stops before opening results')
                if path.startswith('report/'):return self.file(job.folder/'report.html','text/html; charset=utf-8',report=True)
                if path.startswith('excel/'):return self.file(job.folder/'PlasmidTrace_results.xlsx','application/vnd.openxmlformats-officedocument.spreadsheetml.sheet','PlasmidTrace_results.xlsx')
                with state.zip_lock:archive=make_archive(job.folder)
                return self.file(archive,'application/zip','PlasmidTrace_results.zip')
            self.send({'error':'Not found'},404)
        except PermissionError as exc:self.send({'error':str(exc)},403)
        except (KeyError,FileNotFoundError):self.send({'error':'This run is unavailable. Start a new run; previous output folders remain on disk.'},404)
        except (BrokenPipeError,ConnectionResetError):pass
        except (ValueError,OSError) as exc:self.send({'error':str(exc)},400)
    def do_POST(self):
        try:
            path=self.route();origin=self.headers.get('Origin')
            if origin and origin!=f'http://127.0.0.1:{self.server.server_port}':raise PermissionError('Invalid origin')
            length=int(self.headers.get('Content-Length','0'))
            if not 0<=length<=1024*1024:return self.send({'error':'Request too large'},413)
            data=json.loads(self.rfile.read(length) or b'{}');state=self.server.state
            if not isinstance(data,dict):raise ValueError('Request must be a JSON object')
            if path in ('api/setup/detect','api/setup/install'):return self.send({'id':state.setup(data,install=path.endswith('/install'))},202)
            if path=='api/setup/cancel':
                if state.setup_task:threading.Thread(target=state.setup_task.cancel,daemon=True).start()
                return self.send({'ok':True})
            if path=='api/run':return self.send({'id':state.start(data)},202)
            if path=='api/check':return self.send(state.check(data))
            if path=='api/discover':return self.send(state.discover(data))
            if path=='api/pick':return self.send(state.pick(data.get('kind','file')))
            if path.startswith('api/cancel/'):
                job=state.jobs[path.split('/')[-1]];threading.Thread(target=job.cancel,daemon=True).start();return self.send({'ok':True})
            if path.startswith('api/open/'):
                p=state.jobs[path.split('/')[-1]].folder
                if os.name=='nt':os.startfile(p)
                else:subprocess.Popen(['xdg-open',str(p)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
                return self.send({'ok':True})
            self.send({'error':'Not found'},404)
        except PermissionError as exc:self.send({'error':str(exc)},403)
        except KeyError:self.send({'error':'Run not found'},404)
        except (ValueError,OSError,subprocess.SubprocessError) as exc:self.send({'error':str(exc)},400)

def serve(port=0,root=ROOT):
    server=ThreadingHTTPServer(('127.0.0.1',port),Handler);server.token=secrets.token_urlsafe(24);server.state=State(root);return server

def picker(kind):
    import tkinter as tk
    from tkinter import filedialog
    root=tk.Tk();root.withdraw();root.attributes('-topmost',True)
    try:
        path=filedialog.askdirectory(parent=root,initialdir=str(ROOT),title='Choose result or database folder') if kind=='folder' else filedialog.askopenfilename(parent=root,initialdir=str(ROOT),title='Select a file for PlasmidTrace')
        print(json.dumps(path))
    finally:root.destroy()

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--pick',choices=['file','folder']);parser.add_argument('--no-browser',action='store_true');parser.add_argument('--port',type=int,default=0);args=parser.parse_args()
    if args.pick:return picker(args.pick)
    server=serve(args.port);url=f'http://127.0.0.1:{server.server_port}/{server.token}/'
    print('PlasmidTrace · browser app preview 0.3',flush=True);print(url,flush=True);print('Keep this launcher open while analysis runs. Press Ctrl+C to stop.',flush=True)
    if not args.no_browser:threading.Timer(.5,webbrowser.open,args=(url,)).start()
    try:server.serve_forever()
    except KeyboardInterrupt:
        if server.state.setup_task and not server.state.setup_task.finished.is_set():
            server.state.setup_task.cancel();server.state.setup_task.finished.wait(30)
        active=server.state.active
        if active:
            job=server.state.jobs[active]
            if not job.finished.is_set():job.cancel();job.finished.wait(30)
    finally:server.server_close()

if __name__=='__main__':main()
