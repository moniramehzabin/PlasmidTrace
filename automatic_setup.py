"""Coordinate automatic Linux discovery and one-click isolated setup from Windows."""
import json,os,subprocess,threading,time,uuid
from pathlib import Path
from local_app import Backend,ROOT

def assert_f_storage(distro):
    if os.name!='nt':return
    import winreg
    with winreg.OpenKey(winreg.HKEY_CURRENT_USER,r'Software\Microsoft\Windows\CurrentVersion\Lxss') as key:
        for i in range(winreg.QueryInfoKey(key)[0]):
            with winreg.OpenKey(key,winreg.EnumKey(key,i)) as child:
                try:name=winreg.QueryValueEx(child,'DistributionName')[0];base=winreg.QueryValueEx(child,'BasePath')[0]
                except OSError:continue
                if name==distro:
                    if base.removeprefix('\\\\?\\').lower().startswith('f:\\'):return
                    raise ValueError('Automatic downloads require a WSL distribution stored on F:. Selected '+distro+' is stored elsewhere. No files were moved or downloaded.')
    raise ValueError('Selected WSL distribution was not found. Choose your existing Ubuntu installation.')

class SetupTask:
    def __init__(self,distro,root,install=False):
        self.id=uuid.uuid4().hex;self.folder=Path(root)/'setup_runs'/self.id;self.folder.mkdir(parents=True)
        self.backend=Backend('wsl' if os.name=='nt' else 'local',distro,'');self.state='starting';self.process=None;self.finished=threading.Event();self.cancelled=False;self.install=install;self.lock=threading.RLock()
    def run(self):
        try:
            if self.install:assert_f_storage(self.backend.distro)
            result=self.folder/'result.json'
            args=[self.backend.path(ROOT/'setup_worker.py'),'--install' if self.install else '--scan','--result',self.backend.path(result)]
            command=self.backend.wrap(args,self.folder)
            with (self.folder/'console.log').open('wb') as f:
                with self.lock:
                    if self.cancelled:self.state='cancelled';return
                    kw={'start_new_session':True} if os.name!='nt' else {'creationflags':subprocess.CREATE_NEW_PROCESS_GROUP}
                    self.process=subprocess.Popen(command,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,**kw);self.state='installing' if self.install else 'checking'
                code=self.process.wait()
            self.state='cancelled' if self.cancelled else 'finished' if code==0 else 'failed'
        except Exception as exc:
            self.state='failed';(self.folder/'console.log').write_text(str(exc),encoding='utf-8')
        finally:self.finished.set()
    def cancel(self):
        from local_app import Job
        Job.cancel(self)
    def status(self):
        from browser_app import read_tail
        result={}
        try:result=json.loads((self.folder/'result.json').read_text(encoding='utf-8'))
        except (OSError,ValueError):pass
        return dict(id=self.id,state=self.state,done=self.finished.is_set(),log=read_tail(self.folder/'console.log'),result=result,install=self.install)

def apply_runtime(settings,result):
    s=dict(settings)
    if not result.get('ready'):raise ValueError('Setup is not ready. Use Prepare tools and databases first.')
    s.update(backend='wsl' if os.name=='nt' else 'local',environment=result['environment'] if os.name=='nt' else '',
        plasmidfinder_python=result.get('plasmidfinder_python',''),plasmidfinder_executable=result.get('plasmidfinder_executable','plasmidfinder.py'),
        plasmidfinder_db=result['plasmidfinder_db'],mob_db=result['mob_db'],mob_executable=result['tools']['mob_recon'])
    s['runtime_python']=result['python']
    for flag,report in [('run_barrnap','rrna_gff'),('run_mob','mob_report'),('run_plasmidfinder','plasmidfinder_results')]:s[flag]=not bool(s.get(report))
    return s

def validate_screening(s):
    if s['workflow']=='cargo':return
    markers=bool(s.get('run_mob') or s.get('mob_report') or s.get('run_plasmidfinder') or s.get('plasmidfinder_results'))
    if not markers and not s.get('exploratory'):
        raise ValueError('No plasmid-marker checks are enabled. Use Automatic screening and Prepare tools and databases, or explicitly choose an exploratory run in Advanced settings.')
