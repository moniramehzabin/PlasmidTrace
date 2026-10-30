"""Linux-side isolated installation and discovery for PlasmidTrace.

Invoked by the local app, never with sudo. Databases are initialized explicitly.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, json, os, platform, shutil, subprocess, sys, time, urllib.request
from pathlib import Path

PROFILE='core-1'
PACKAGES=['python=3.11','mob_suite=3.1.9','unicycler=0.5.0','barrnap','blast','bwa','samtools','mash','kma','git','pip']
PF_PACKAGE='PlasmidFinder==3.0.3'
CGE_PACKAGE='cgecore==1.5.6'
PF_DATABASE='https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git'
TOOLS=['barrnap','mob_recon','mob_init','unicycler','spades','blastn','makeblastdb','bwa','samtools','mash','kma_index']
MOB_FILES=['status.txt','ncbi_plasmid_full_seqs.fas','ncbi_plasmid_full_seqs.fas.msh','clusters.txt','taxa.sqlite']

def default_runtime():return Path.home()/'.local/share/plasmidtrace/runtime-core-1'
def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def write_json(path,data):
    path=Path(path);temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(data,indent=2),encoding='utf-8');temp.replace(path)
def run(argv,env=None,cwd=None,timeout=None):
    print('Running: '+' '.join(map(str,argv)),flush=True)
    subprocess.run(list(map(str,argv)),env=env,cwd=cwd,check=True,timeout=timeout)
def conda_executable():
    for root in [Path.home()/'miniforge3',Path.home()/'miniconda3',Path.home()/'anaconda3',Path('/opt/conda')]:
        if (root/'bin/conda').is_file():return root/'bin/conda'
    found=shutil.which('conda')
    return Path(found) if found else None

def inspect(runtime):
    runtime=Path(runtime);tools={name:shutil.which(name) for name in TOOLS}
    missing=[name for name,path in tools.items() if not path]
    pf_module=importlib.util.find_spec('plasmidfinder') is not None
    pf_exe=shutil.which('plasmidfinder.py')
    if not pf_module and not pf_exe:missing.append('PlasmidFinder')
    pf_paths=[runtime/'databases/plasmidfinder',Path.home()/'databases/plasmidfinder_db',Path.home()/'plasmidfinder_db',Path('/opt/plasmidfinder_db')]
    mob_paths=[runtime/'databases/mob']
    spec=importlib.util.find_spec('mob_suite')
    if spec and spec.submodule_search_locations:
        mob_paths += [Path(p)/'databases' for p in spec.submodule_search_locations]
    pf_db=next((p for p in pf_paths if (p/'config').is_file() and any(f.stat().st_size for f in p.glob('*.fsa'))),None)
    mob_db=next((p for p in mob_paths if all((p/n).is_file() and (p/n).stat().st_size>0 for n in MOB_FILES) and any(p.glob('*.nhr'))),None)
    if not pf_db:missing.append('PlasmidFinder database')
    if not mob_db:missing.append('MOB-suite database')
    if pf_db==runtime/'databases/plasmidfinder' and not (runtime/'checks/positive_control.json').is_file():missing.append('Database positive-control verification')
    return dict(ready=not missing,profile=PROFILE,missing=missing,tools=tools,python=sys.executable,
        environment=str(Path(sys.executable).resolve().parent.parent),plasmidfinder_python=sys.executable if pf_module else '',
        plasmidfinder_executable=pf_exe or 'plasmidfinder.py',plasmidfinder_db=str(pf_db or ''),mob_db=str(mob_db or ''),runtime=str(runtime))

def scan(runtime,worker=None):
    """Probe system and existing environments automatically, with bounded waits."""
    runtime=Path(runtime);worker=Path(worker or __file__).resolve()
    candidates=[runtime/'env/bin/python',Path(sys.executable)]
    conda=conda_executable()
    if conda:
        try:
            data=json.loads(subprocess.check_output([str(conda),'env','list','--json'],text=True,timeout=30))
            candidates += [Path(p)/'bin/python' for p in data.get('envs',[])[:16]]
        except (OSError,ValueError,subprocess.SubprocessError):pass
    seen=set();best=None
    for executable in candidates:
        if not executable.is_file() or str(executable) in seen:continue
        seen.add(str(executable));print('Checking installed environment: '+str(executable.parent.parent),flush=True)
        try:
            env=dict(os.environ,PATH=str(executable.parent)+os.pathsep+os.environ.get('PATH',''))
            output=subprocess.check_output([str(executable),str(worker),'--inspect','--runtime',str(runtime)],text=True,env=env,stderr=subprocess.DEVNULL,timeout=25)
            result=json.loads(output)
            if best is None or len(result['missing'])<len(best['missing']):best=result
            if result['ready']:return result
        except (OSError,ValueError,subprocess.SubprocessError):continue
    return best or dict(ready=False,missing=['Linux Python and analysis tools'],runtime=str(runtime),profile=PROFILE)

def bootstrap_conda(runtime,env):
    found=conda_executable()
    if found:return found
    machine=platform.machine()
    if machine not in ('x86_64','aarch64'):raise ValueError('Automatic installation supports Linux x86_64 and aarch64.')
    filename=f'Miniforge3-Linux-{machine}.sh';base='https://github.com/conda-forge/miniforge/releases/latest/download/'
    installer=runtime/'downloads'/filename;installer.parent.mkdir(parents=True,exist_ok=True)
    print('Downloading Miniforge installer and its published SHA-256 checksum…',flush=True)
    urllib.request.urlretrieve(base+filename,installer)
    with urllib.request.urlopen(base+filename+'.sha256',timeout=60) as response:checksum=response.read().decode().split()[0]
    if digest(installer)!=checksum:raise ValueError('Installer checksum mismatch. Setup stopped.')
    prefix=runtime/'miniforge'
    if prefix.exists():prefix.rename(runtime/('miniforge.incomplete.'+str(time.time_ns())))
    run(['bash',installer,'-b','-p',prefix],env=env)
    return prefix/'bin/conda'

def first_record(path):
    chunks=[];name=None
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                if name:break
                name=line[1:].strip().split()[0]
            elif name:chunks.append(line.strip())
    seq=''.join(chunks)
    if not name or not seq:raise ValueError('Database FASTA has no usable sequence: '+str(path))
    return name,seq

def install(runtime):
    if platform.system()!='Linux':raise ValueError('Setup runs inside Linux / WSL, not Windows Python.')
    runtime=Path(runtime).expanduser().resolve()
    if runtime==Path('/') or runtime==Path.home():raise ValueError('Choose the dedicated app runtime directory.')
    runtime.mkdir(parents=True,exist_ok=True)
    if shutil.disk_usage(runtime).free < 8*1024**3:raise ValueError('At least 8 GB of free space is required for first-time tools and databases.')
    if hasattr(os,'sched_getaffinity'):os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:2]))
    import fcntl
    with (runtime/'setup.lock').open('w') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('Another setup is already using this runtime. Let it finish first.')
        env=dict(os.environ)
        for key,relative in [('CONDA_PKGS_DIRS','cache/conda'),('PIP_CACHE_DIR','cache/pip'),('TMPDIR','tmp'),('XDG_CACHE_HOME','cache'),('TMP','tmp'),('TEMP','tmp')]:
            p=runtime/relative;p.mkdir(parents=True,exist_ok=True);env[key]=str(p)
        conda=bootstrap_conda(runtime,env);prefix=runtime/'env'
        print('STEP 1/4 · Installing isolated analysis tools. Existing environments will not be changed.',flush=True)
        operation='install' if (prefix/'conda-meta').is_dir() else 'create'
        if operation=='create' and prefix.exists():prefix.rename(runtime/('env.incomplete.'+str(time.time_ns())))
        run([conda,operation,'--yes','--prefix',prefix,'--override-channels','-c','conda-forge','-c','bioconda','--strict-channel-priority']+PACKAGES,env=env)
        python=prefix/'bin/python';env['PATH']=str(prefix/'bin')+os.pathsep+env.get('PATH','')
        constraints=runtime/'pip-constraints.txt'
        constraints.write_text('numpy<1.23.5\npandas<=1.5.3\nbiopython<2\ncgecore==1.5.6\n',encoding='utf-8')
        run([python,'-m','pip','install','--constraint',constraints,PF_PACKAGE,CGE_PACKAGE],env=env)
        run([python,'-m','pip','check'],env=env)
        if not (prefix/'bin/spades').exists() and (prefix/'bin/spades.py').exists():(prefix/'bin/spades').symlink_to('spades.py')
        print('STEP 2/4 · Downloading and indexing the PlasmidFinder database.',flush=True)
        db=runtime/'databases/plasmidfinder';db.parent.mkdir(parents=True,exist_ok=True)
        if not (db/'config').is_file():
            stage=runtime/'databases/plasmidfinder.partial'
            if stage.exists():stage.rename(stage.with_name('plasmidfinder.partial.previous.'+str(time.time_ns())))
            run([prefix/'bin/git','clone','--depth','1',PF_DATABASE,stage],env=env)
            stage.rename(db)
        if not list(db.glob('*.fsa')):raise ValueError('PlasmidFinder database is empty.')
        installer=db/'INSTALL.py'
        if installer.is_file():run([python,installer,prefix/'bin/kma_index'],env=env,cwd=db)
        for fasta in sorted(db.glob('*.fsa')):
            if not Path(str(fasta)+'.nhr').exists():run([prefix/'bin/makeblastdb','-in',fasta,'-dbtype','nucl'],env=env)
        print('STEP 3/4 · Initializing MOB-suite databases. This is a large first-time download.',flush=True)
        mob=runtime/'databases/mob';mob.mkdir(exist_ok=True)
        if not all((mob/name).is_file() and (mob/name).stat().st_size>0 for name in MOB_FILES):
            # This directory is owned exclusively by this profile; our setup lock
            # excludes concurrent initializers. A cancelled mob_init can leave this lock.
            (mob/'.lock').unlink(missing_ok=True)
            run([prefix/'bin/mob_init','-d',mob],env=env,cwd=runtime)
        print('STEP 4/4 · Checking database files, installed programs and PlasmidFinder positive control.',flush=True)
        control=runtime/'checks'/str(time.time_ns());control.mkdir(parents=True)
        fasta=next(p for p in sorted(db.glob('*.fsa')) if p.stat().st_size)
        _,sequence=first_record(fasta);query=control/'marker.fasta';query.write_text('>setup_control\n'+sequence+'\n')
        output=control/'plasmidfinder';output.mkdir()
        run([python,'-m','plasmidfinder','-i',query,'-o',output,'-p',db,'-x','-t','0.95','-l','0.60'],env=env)
        from finder_adapter import import_hits
        hits=import_hits(output/'results_tab.tsv',{'setup_control':sequence},.95,.60)
        if not hits['setup_control']:raise ValueError('PlasmidFinder positive control failed. Setup is not marked ready.')
        write_json(runtime/'checks/positive_control.json',{'status':'passed','version':'3.0.3','checked_at':time.time(),'query_sha256':digest(query)})
        result=json.loads(subprocess.check_output([str(python),str(Path(__file__).resolve()),'--inspect','--runtime',str(runtime)],env=env,text=True))
        if not result['ready']:raise ValueError('Setup verification failed: '+', '.join(result['missing']))
        result['created_at']=time.time();result['plasmidfinder_version']='3.0.3';result['positive_control']='passed'
        result['plasmidfinder_database_commit']=subprocess.check_output([str(prefix/'bin/git'),'-C',str(db),'rev-parse','HEAD'],env=env,text=True).strip()
        result['mob_database_status_sha256']=digest(mob/'status.txt')
        result['packages']=json.loads(subprocess.check_output([str(conda),'list','--prefix',str(prefix),'--json'],env=env,text=True))
        write_json(runtime/'ready.json',result)
        print('Setup verified. Marker checks can now run automatically.',flush=True)
        return result

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--runtime',type=Path,default=default_runtime());parser.add_argument('--inspect',action='store_true');parser.add_argument('--scan',action='store_true');parser.add_argument('--install',action='store_true');parser.add_argument('--result',type=Path)
    a=parser.parse_args()
    try:
        data=inspect(a.runtime) if a.inspect else scan(a.runtime) if a.scan else install(a.runtime) if a.install else None
        if data is None:parser.error('Choose an operation')
        if a.result:write_json(a.result,data)
        if a.inspect:print(json.dumps(data))
        return 0
    except Exception as exc:
        print('SETUP ERROR: '+str(exc),flush=True)
        if a.result:write_json(a.result,{'ready':False,'error':str(exc)})
        return 1
if __name__=='__main__':raise SystemExit(main())
