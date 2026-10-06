import unittest,tempfile,os,subprocess,sys,json
from pathlib import Path
class OrchestrationTest(unittest.TestCase):
    def test_active_workflow_with_simulated_executables(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);binpath=root/'bin';binpath.mkdir()
            fixture=root/'input.fa';fixture.write_text('>12\n'+'A'*1000+'\n')
            graph=root/'input.gfa';graph.write_text('S\t12\t'+'A'*1000+'\nL\t12\t+\t12\t+\t0M\n')
            for name in ['unicycler','spades','bwa','samtools','blastn','makeblastdb']:
                exe=binpath/name
                exe.write_text('''#!/usr/bin/env python3
import sys,os
from pathlib import Path
name=Path(sys.argv[0]).name;a=sys.argv[1:]
if '--version' in a or '-version' in a or not a: print('simulated');sys.exit(0)
if name=='unicycler':
    p=Path(a[a.index('-o')+1]);p.mkdir()
    (p/'assembly.fasta').write_text('>12\\n'+'A'*1000+'\\n')
    (p/'assembly.gfa').write_text('S\\t12\\t'+'A'*1000+'\\nL\\t12\\t+\\t12\\t+\\t0M\\n')
elif name=='blastn':
    Path(a[a.index('-out')+1]).write_text('12\\t12\\t100\\t1\\t1000\\t1000\\t1000\\t2000\\n')
elif name=='bwa' and a[0]=='mem':print('simulated mapping stream')
elif name=='samtools' and a[0]=='sort':
    sys.stdin.read();Path(a[a.index('-o')+1]).write_bytes(b'simulated')
elif name=='samtools' and a[0]=='view':
    region=a[-1];ref=region.split(':')[0]
    print('read1\\t0\\t'+ref+'\\t451\\t30\\t100M\\t*\\t0\\t0\\t'+'A'*100+'\\t'+'I'*100+'\\tNM:i:0')
''');exe.chmod(0o755)
            r1=root/'r1';r2=root/'r2';r1.write_text('fixture');r2.write_text('fixture')
            env=os.environ.copy();env['PATH']=str(binpath)+os.pathsep+env['PATH']
            script=Path(__file__).with_name('reconstruct.py')
            r=subprocess.run([sys.executable,str(script),'--assembly',str(fixture),'--gfa',str(graph),'--reads1',str(r1),'--reads2',str(r2),'--out',str(root/'out')],env=env,capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            data=json.loads((root/'out/workflow.json').read_text())
            self.assertEqual(data['state'],'completed');self.assertEqual(data['results'][0]['crossing_reads'],1)
            self.assertEqual(data['results'][0]['mapping_reference'],'junction only')
            self.assertFalse(list((root/'out').glob('junction_*/competitive_reference.fasta')))
            self.assertTrue(any('unicycler' in cmd[0] for cmd in data['commands']))
if __name__=='__main__':unittest.main()
