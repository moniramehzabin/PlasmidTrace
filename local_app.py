#!/usr/bin/env python3
"""PlasmidTrace local analysis runner. Standard library only; engine unchanged."""
from __future__ import annotations
import json, os, re, shlex, shutil, signal, subprocess, sys, threading, time, uuid, webbrowser, zipfile
from pathlib import Path, PureWindowsPath
import html
ROOT = Path(__file__).resolve().parent
APP_VERSION = '0.4.0-preview'
CONDA_INIT = '''for p in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3"; do
 if [ -f "$p/etc/profile.d/conda.sh" ]; then . "$p/etc/profile.d/conda.sh"; break; fi
done
'''

class Backend:
    def __init__(self, kind='local', distro='Ubuntu-24.04', environment='', python=None):
        if kind not in ('local','wsl'): raise ValueError('Invalid backend')
        self.kind, self.distro, self.environment = kind, distro.strip(), environment.strip()
        self.python=python or sys.executable
        self.custom_python=bool(python)
        if kind == 'wsl' and not self.distro: raise ValueError('Select a WSL distribution')
    def path(self, value):
        value = str(value).strip()
        if self.kind == 'wsl' and re.match(r'^[A-Za-z]:[\\/]', value):
            p = PureWindowsPath(value)
            return '/mnt/' + p.drive[0].lower() + '/' + '/'.join(p.parts[1:])
        if self.kind == 'wsl' and value.startswith('\\\\'):
            raise ValueError('Use a drive path (such as E:) or a Linux path, not a network path')
        return value
    def shell(self, script):
        if os.name == 'nt':
            # Pass a script file, avoiding Windows/WSL command-line quote loss.
            folder = ROOT / 'setup_runs' / 'shell'
            folder.mkdir(parents=True, exist_ok=True)
            path = folder / (uuid.uuid4().hex + '.sh')
            path.write_bytes(script.replace('\r\n', '\n').encode('utf-8'))
            return ['wsl.exe', '-d', self.distro, '--exec', 'bash', self.path(path)]
        return ['wsl.exe','-d',self.distro,'--','bash','-lc',script]
    def prefix(self):
        if self.kind=='wsl' and self.custom_python and self.python.startswith('/'):
            return 'export PATH='+shlex.quote(str(Path(self.python).parent))+':"$PATH"\nPT_PYTHON='+shlex.quote(self.python)+'\n'
        return CONDA_INIT + ('conda activate '+shlex.quote(self.environment)+' || exit $?\n' if self.environment else '') + 'PT_PYTHON=$(command -v python3 || command -v python) || exit 127\n'
    def wrap(self, arguments, folder):
        if self.kind == 'local': return [self.python, '-u'] + list(arguments)
        pid = self.path(folder / 'linux_process_group.pid')
        # setsid gives this run its own Linux process group for cancellation.
        marker = self.path(folder / 'cancel.requested')
        guard = '[ ! -f '+shlex.quote(marker)+' ] || exit 130\n'
        script = self.prefix() + guard + 'setsid "$PT_PYTHON" ' + shlex.join(['-u']+list(arguments)) + ' &\npid=$!\nprintf "%s" "$pid" > '+shlex.quote(pid)+'\nif [ -f '+shlex.quote(marker)+' ]; then kill -TERM -- "-$pid" 2>/dev/null; fi\nwait "$pid"\n'
        return self.shell(script)
    def probe(self, names=None):
        code = 'import json,shutil,sys;print(json.dumps({"python":sys.executable,"tools":{n:shutil.which(n) for n in '+repr(names or ['barrnap','mob_recon','plasmidfinder.py','unicycler','spades','blastn','makeblastdb','bwa','samtools','amrfinder'])+'}}))'
        cmd = [self.python,'-c',code] if self.kind=='local' else self.shell(self.prefix()+'"$PT_PYTHON" -c '+shlex.quote(code))
        result = subprocess.run(cmd,capture_output=True,text=True,timeout=45,env=dict(os.environ,PATH=str(Path(self.python).parent)+os.pathsep+os.environ.get('PATH','')) if self.kind=='local' else None)
        if result.returncode: raise ValueError((result.stderr or result.stdout).strip())
        return json.loads(result.stdout.strip().splitlines()[-1])

PATH_FLAGS = {'assembly':'--assembly','rrna_gff':'--rrna-gff','mob_report':'--mob-report','gfa':'--gfa','annotation_genbank':'--annotation-genbank','annotation_gff':'--annotation-gff','annotation_fasta':'--annotation-fasta','plasmidfinder_results':'--plasmidfinder-results','plasmidfinder_db':'--plasmidfinder-db','mob_db':'--mob-db'}

def number(settings, key, minimum, maximum):
    try: value=int(settings.get(key,minimum))
    except (TypeError,ValueError): raise ValueError(key+' must be a whole number')
    if not minimum<=value<=maximum: raise ValueError(f'{key} must be between {minimum} and {maximum}')
    return str(value)

def cargo_config(s, backend):
    config={}
    for tool in ('amrfinder','virulencefinder','geneminer'):
        mode=s.get(tool+'_mode','off')
        if mode=='off': continue
        if mode not in ('import','run'): raise ValueError('Invalid cargo mode')
        def required(name):
            v=s.get(tool+'_'+name,'').strip()
            if not v: raise ValueError(tool+': choose '+name.replace('_',' '))
            return backend.path(v)
        if mode=='import':
            item={'result':required('result')}
        elif tool=='amrfinder':
            item={'command':[backend.path(s.get('amrfinder_exe','amrfinder')),'-n','{assembly}','--plus','--threads',number(s,'threads',1,256),'-o','{out}/amrfinder.tsv'],'result':'{out}/amrfinder.tsv'}
            if s.get('amrfinder_db'): item['command']+=['--database',backend.path(s['amrfinder_db'])]
        elif tool=='virulencefinder':
            item={'command':[backend.path(s.get('virulencefinder_python','python')),'-m','virulencefinder','-ifa','{assembly}','-o','{out}','-p',required('db'),'-b',backend.path(s.get('virulencefinder_blast','blastn')),'-d','all','-t','0.90','-l','0.60','-x'],'result':'{out}/results_tab.tsv'}
        else:
            item={'command':[backend.path(s.get('geneminer_python','python')),required('script'),required('genbank'),'--reference-db',required('db'),'--reference-metadata',required('metadata'),'--annotation-rules',required('rules'),'--threads',number(s,'threads',1,256),'--outdir','{out}'],'result':'{out}/Bioremediation_Gene_Miner_Report_all_candidates.tsv','workbook':'{out}/Bioremediation_Gene_Miner_Report.xlsx','provenance_files':[required('script'),required('db'),required('metadata'),required('rules')]}
        if tool=='geneminer':
            item['genbank']=required('genbank')
            if mode=='import' and s.get('geneminer_workbook'): item['workbook']=backend.path(s['geneminer_workbook'])
        config[tool]=item
    return config

def build_command(s, backend, folder):
    mode=s.get('workflow','screen')
    if mode not in ('screen','reconstruct','reuse','recover','cargo'): raise ValueError('Choose a workflow')
    if not s.get('assembly','').strip(): raise ValueError('Choose an assembly FASTA')
    if s.get('run_barrnap') and s.get('rrna_gff'): raise ValueError('Choose barrnap OR an existing rRNA GFF')
    if s.get('run_mob') and s.get('mob_report'): raise ValueError('Choose MOB-Recon OR an existing MOB report')
    if s.get('run_plasmidfinder') and s.get('plasmidfinder_results'): raise ValueError('Choose PlasmidFinder OR an existing results table')
    if s.get('annotation_genbank') and s.get('annotation_gff'): raise ValueError('Choose GenBank OR GFF annotation')
    if s.get('annotation_fasta') and not s.get('annotation_gff'): raise ValueError('Annotation FASTA requires annotation GFF')
    script='plasmid_evidence.py' if mode=='screen' else 'cargo_analysis.py' if mode=='cargo' else 'reconstruct.py'
    args=[backend.path(ROOT/script),'--assembly',backend.path(s['assembly']),'--out',backend.path(folder/'output')]
    if mode!='cargo':
        args+=['--threads',number(s,'threads',1,256)]
        for key,flag in PATH_FLAGS.items():
            if key!='assembly' and s.get(key): args += [flag,backend.path(s[key])]
        for key in ('run_barrnap','run_mob','run_plasmidfinder'):
            if s.get(key): args+=['--'+key.replace('_','-')]
        if s.get('run_plasmidfinder') and not s.get('plasmidfinder_db'): raise ValueError('Choose the PlasmidFinder database folder')
        for key in ('mob_executable','plasmidfinder_executable','plasmidfinder_python'):
            if s.get(key): args+=['--'+key.replace('_','-'),backend.path(s[key])]
    if mode in ('reconstruct','reuse','recover'):
        for key in ('reads1','reads2'):
            if not s.get(key): raise ValueError('Choose both paired read files')
            args+=['--'+key,backend.path(s[key])]
        if backend.path(s['reads1'])==backend.path(s['reads2']): raise ValueError('Read mates must be different files')
        args+=['--memory-gb',number(s,'memory',1,1024)]
        if mode in ('reuse','recover'):
            if not s.get('previous'): raise ValueError('Choose the previous run folder')
            args += ['--reuse-unicycler' if mode=='reuse' else '--recover-from',backend.path(s['previous'])]
        if mode=='recover': args+=['--unicycler-python',backend.path(s.get('recovery_python') or 'python')]
    config=cargo_config(s,backend)
    if mode=='screen' and config: raise ValueError('Use reconstruction or the candidate gene-analysis workflow for cargo modules')
    if mode=='cargo' and not config: raise ValueError('Enable at least one gene-analysis module')
    if config:
        p=folder/'cargo_config.json';p.write_text(json.dumps(config,indent=2),encoding='utf-8')
        args+=['--config' if mode=='cargo' else '--cargo-config',backend.path(p)]
    return backend.wrap(args,folder)

def required_tools(s):
    names=[]
    if s.get('workflow')!='cargo':
        if s.get('run_barrnap'): names+=['barrnap']
        if s.get('run_mob'): names += [s.get('mob_executable') or 'mob_recon']
        if s.get('run_plasmidfinder'): names += [s.get('plasmidfinder_python') or s.get('plasmidfinder_executable') or 'plasmidfinder.py']
    if s.get('workflow') in ('reconstruct','reuse','recover'):
        names+=['bwa','samtools','blastn','makeblastdb']
        if s['workflow']!='reuse': names+=['unicycler','spades']
    for tool,key,default in [('amrfinder','amrfinder_exe','amrfinder'),('virulencefinder','virulencefinder_python','python'),('geneminer','geneminer_python','python')]:
        if s.get(tool+'_mode')=='run':names.append(s.get(key) or default)
    return list(dict.fromkeys(names))

def export_report(folder, state):
    from app_export import export_results
    from assessment import summarize
    try: settings=json.loads((folder/'app_run.json').read_text())['settings']
    except (OSError,ValueError,KeyError):settings={}
    assessment=summarize(folder,settings)
    (folder/'assessment.json').write_text(json.dumps(assessment,indent=2),encoding='utf-8')
    sections=['<h2>'+html.escape(assessment['label'])+'</h2><ul>'+''.join('<li>'+html.escape(w)+'</li>' for w in assessment['warnings'])+'</ul>']
    for p in sorted((folder/'output').rglob('report.md')):
        sections.append('<h2>'+html.escape(str(p.relative_to(folder/'output')))+'</h2><pre>'+html.escape(p.read_text(errors='replace'))+'</pre>')
    (folder/'report.html').write_text('<!doctype html><meta charset="utf-8"><title>PlasmidTrace results</title><style>body{font:16px/1.6 system-ui;max-width:1000px;margin:40px auto;padding:20px;color:#173c39}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#eff5f3;padding:20px;border-radius:12px}</style><h1>PlasmidTrace</h1><p>Run status: '+html.escape(state)+'</p><p>Computational evidence for review. Candidate calls and read support do not confirm plasmid identity. Predicted genes do not establish activity.</p>'+''.join(sections),encoding='utf-8')
    export_results(folder)

def make_archive(folder):
    archive=folder/'PlasmidTrace_results.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(folder.rglob('*')):
            if p.is_file() and not p.is_symlink() and p!=archive: z.write(p,p.relative_to(folder))
    return archive

class Job:
    def __init__(self,s,backend,root):
        self.settings=dict(s);self.backend=backend
        self.folder=Path(root).expanduser().resolve()/('run_'+time.strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:8])
        self.folder.mkdir(parents=True)
        self.command=build_command(s,backend,self.folder)
        self.process=None;self.state='starting';self.cancelled=False;self.lock=threading.RLock();self.finished=threading.Event()
        (self.folder/'app_run.json').write_text(json.dumps({'app_version':APP_VERSION,'engine_version':'1.0.0','settings':s,'backend':{'kind':backend.kind,'distro':backend.distro,'environment':backend.environment},'command':self.command,'state':self.state},indent=2),encoding='utf-8')
    def run(self):
        try:
            with (self.folder/'console.log').open('wb') as log:
                names = required_tools(self.settings)
                if names:
                    probe = self.backend.probe([self.backend.path(n) for n in names])
                    missing = [k for k,v in probe['tools'].items() if not v]
                    if missing: raise ValueError('Missing tools in the selected environment: '+', '.join(missing)+'. Expand the run environment in Run settings and choose the correct installed environment.')
                with self.lock:
                    if self.cancelled: self.state='cancelled';return
                    kw={'start_new_session':True} if os.name!='nt' else {'creationflags':subprocess.CREATE_NEW_PROCESS_GROUP}
                    self.process=subprocess.Popen(self.command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,env=dict(os.environ,PYTHONUNBUFFERED='1',PATH=str(Path(self.backend.python).parent)+os.pathsep+os.environ.get('PATH','')) if self.backend.kind=='local' else dict(os.environ,PYTHONUNBUFFERED='1'),**kw)
                    self.state='running'
                code=self.process.wait()
                final_state='cancelled' if self.cancelled else 'finished' if code==0 else 'failed'
                self.state='exporting'
            try: export_report(self.folder,final_state)
            except Exception as exc:
                with (self.folder/'console.log').open('a') as f:f.write('\nResult export warning: '+str(exc))
            self.state=final_state
        except Exception as exc:
            self.state='cancelled' if self.cancelled else 'failed'
            (self.folder/'console.log').write_text(str(exc),encoding='utf-8')
        finally:
            p=self.folder/'app_run.json';data=json.loads(p.read_text());data.update(state=self.state,exit_code=self.process.returncode if self.process else None);p.write_text(json.dumps(data,indent=2),encoding='utf-8');self.finished.set()
    def cancel(self):
        with self.lock:
            if self.state not in ('starting','running','checking','installing'): return
            self.cancelled=True
            (self.folder/'cancel.requested').touch()
            if self.backend.kind=='wsl':
                # Never kill an entire WSL distribution. Stop only this run's group.
                pid=self.backend.path(self.folder/'linux_process_group.pid')
                script='if [ -f '+shlex.quote(pid)+' ]; then p=$(cat '+shlex.quote(pid)+'); case "$p" in ""|*[!0-9]*) exit 1;; esac; kill -TERM -- "-$p" 2>/dev/null; sleep 2; kill -KILL -- "-$p" 2>/dev/null; fi'
                subprocess.run(self.backend.shell(script),capture_output=True,timeout=15)
            if self.process and self.process.poll() is None:
                if os.name=='nt':subprocess.run(['taskkill','/PID',str(self.process.pid),'/T','/F'],capture_output=True)
                else:
                    try:os.killpg(self.process.pid,signal.SIGTERM)
                    except ProcessLookupError:pass
                    def force():
                        time.sleep(2)
                        try:os.killpg(self.process.pid,signal.SIGKILL)
                        except ProcessLookupError:pass
                    threading.Thread(target=force,daemon=True).start()

# Compatibility entry point: the browser app replaces the desktop window.
if __name__ == '__main__':
    from browser_app import main
    main()
