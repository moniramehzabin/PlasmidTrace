import json,os,sys,tempfile,threading,unittest
from pathlib import Path
from unittest.mock import patch
from automatic_setup import SetupTask,apply_runtime,validate_screening
from browser_app import State,clean_settings
from local_app import Backend,build_command
from assessment import summarize
import setup_worker

class AutomaticSetupTests(unittest.TestCase):
    def runtime(self):
        return dict(ready=True,environment='/runtime/env',python='/runtime/env/bin/python',tools={'mob_recon':'/runtime/env/bin/mob_recon'},plasmidfinder_python='/runtime/env/bin/python',plasmidfinder_executable='plasmidfinder.py',plasmidfinder_db='/runtime/databases/plasmidfinder',mob_db='/runtime/databases/mob')
    def test_no_silent_empty_screening(self):
        s=clean_settings({'automatic':False})
        with self.assertRaisesRegex(ValueError,'No plasmid-marker'):validate_screening(s)
        s['exploratory']=True;validate_screening(s)
        s['exploratory']=False;s['plasmidfinder_results']='saved.tsv';validate_screening(s)
    def test_default_automatic_mode_requires_ready_runtime(self):
        with tempfile.TemporaryDirectory() as d:
            state=State(d)
            with self.assertRaises(ValueError):state.start({'assembly':'input.fa'})
            self.assertFalse(state.jobs)
    def test_setup_failure_cannot_be_applied(self):
        with self.assertRaisesRegex(ValueError,'not ready'):apply_runtime(clean_settings({}),{'ready':False})
    def test_automatic_preset_enables_all_checks_and_database_paths(self):
        s=apply_runtime(clean_settings({}),self.runtime())
        for key in ('run_barrnap','run_mob','run_plasmidfinder'):self.assertTrue(s[key])
        self.assertEqual(s['mob_db'],'/runtime/databases/mob')
        self.assertEqual(s['runtime_python'],'/runtime/env/bin/python')
    def test_existing_reports_replace_matching_checks(self):
        s=apply_runtime(clean_settings({'rrna_gff':'r.gff','mob_report':'mob.tsv','plasmidfinder_results':'pf.tsv'}),self.runtime())
        for key in ('run_barrnap','run_mob','run_plasmidfinder'):self.assertFalse(s[key])
    def test_custom_wsl_runtime_does_not_need_conda_on_path(self):
        b=Backend('wsl',environment='/env',python='/env/bin/python')
        self.assertNotIn('conda activate',b.prefix());self.assertIn('/env/bin/python',b.prefix())
        with tempfile.TemporaryDirectory() as d:
            args=build_command(dict(assembly='a.fa',workflow='screen',threads='2',mob_db='/data/mob',run_mob=True),b,Path(d))
            self.assertIn('--mob-db /data/mob',args[-1]);self.assertIn('"$PT_PYTHON"',args[-1])
    def test_limited_assessment_flags_original_failure_case(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'output').mkdir()
            (root/'output/evidence.tsv').write_text('contig\tstatus\tcoverage_ratio\trrna_fraction\tmob_classification\tplasmidfinder_markers\nscaffold_1\tno positive evidence in supplied checks\tnot assessed\tnot assessed\tnot assessed\tnot assessed\n')
            r=summarize(root,{})
            self.assertTrue(r['limited']);self.assertEqual(r['review_candidates'],0)
            self.assertIn('No plasmid-marker',r['warnings'][0])
    def test_zero_hits_are_assessed_when_checks_really_ran(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'output').mkdir()
            (root/'output/evidence.tsv').write_text('contig\tstatus\tcoverage_ratio\trrna_fraction\tmob_classification\tplasmidfinder_markers\nc\tno positive evidence in supplied checks\t1.0\t0\tchromosome\t\n')
            r=summarize(root,{})
            self.assertFalse(r['limited']);self.assertTrue(r['checks']['PlasmidFinder']);self.assertEqual(r['review_candidates'],0)
    def test_empty_database_directory_never_ready(self):
        with tempfile.TemporaryDirectory() as d:
            with patch('setup_worker.shutil.which',return_value='/tool'),patch('setup_worker.importlib.util.find_spec',return_value=None):
                r=setup_worker.inspect(Path(d))
            self.assertFalse(r['ready']);self.assertIn('MOB-suite database',r['missing']);self.assertIn('PlasmidFinder database',r['missing'])
    def test_installer_has_fixed_isolated_package_plan(self):
        self.assertIn('python=3.11',setup_worker.PACKAGES);self.assertIn('mob_suite=3.1.9',setup_worker.PACKAGES)
        self.assertEqual(setup_worker.PF_PACKAGE,'PlasmidFinder==3.0.3');self.assertTrue(setup_worker.PF_DATABASE.startswith('https://bitbucket.org/genomicepidemiology/'))
    def test_missing_database_control_cannot_be_ready(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);pf=root/'databases/plasmidfinder';pf.mkdir(parents=True);(pf/'config').write_text('db\tdesc\tref\n');(pf/'db.fsa').write_text('>x\nAAAA\n')
            mob=root/'databases/mob';mob.mkdir()
            for n in setup_worker.MOB_FILES:(mob/n).write_text('test')
            (mob/'x.nhr').write_text('index')
            with patch('setup_worker.shutil.which',return_value='/tool'),patch('setup_worker.importlib.util.find_spec',return_value=None):r=setup_worker.inspect(root)
            self.assertFalse(r['ready']);self.assertIn('Database positive-control verification',r['missing'])
    def test_setup_and_analysis_cannot_run_concurrently(self):
        with tempfile.TemporaryDirectory() as d:
            state=State(d)
            class Running:
                finished=threading.Event()
            state.setup_task=Running()
            with self.assertRaisesRegex(ValueError,'setup finish'):state.start({'automatic':False,'exploratory':True,'assembly':'x.fa'})
    def test_changed_distribution_invalidates_detection(self):
        with tempfile.TemporaryDirectory() as d:
            state=State(d);state.runtime_result=self.runtime();state.setup_distro='Ubuntu-24.04'
            with self.assertRaisesRegex(ValueError,'detection'):state.start({'assembly':'x.fa','distro':'Debian-F'})
    def test_first_database_record_for_positive_control(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'db.fa';p.write_text('>marker1 description\nACGT\nACGT\n>marker2\nTTTT\n')
            self.assertEqual(setup_worker.first_record(p),('marker1','ACGTACGT'))

if __name__=='__main__':unittest.main()

class InstallerOrchestrationTests(unittest.TestCase):
    def exercise(self,root,positive):
        calls=[]
        def fake_run(argv,env=None,cwd=None,timeout=None):
            argv=list(map(str,argv));calls.append((argv,env))
            if 'create' in argv or 'install' in argv and '--prefix' in argv:
                prefix=Path(argv[argv.index('--prefix')+1]);(prefix/'conda-meta').mkdir(parents=True);(prefix/'bin').mkdir()
                for name in setup_worker.TOOLS+['python','git']:(prefix/'bin'/name).touch()
            elif 'clone' in argv:
                p=Path(argv[-1]);p.mkdir();(p/'config').write_text('db\tdescription\tref\n');(p/'db.fsa').write_text('>marker\nACGTACGT\n')
            elif Path(argv[0]).name=='makeblastdb':Path(argv[argv.index('-in')+1]+'.nhr').write_text('index')
            elif Path(argv[0]).name=='mob_init':
                p=Path(argv[-1]);p.mkdir(exist_ok=True)
                for name in setup_worker.MOB_FILES:(p/name).write_text('ready')
                (p/'index.nhr').write_text('index')
            elif '-m' in argv and 'plasmidfinder' in argv:
                p=Path(argv[argv.index('-o')+1]);p.mkdir(exist_ok=True)
                text='Plasmid\tIdentity\tQuery / Template length\tContig\n'
                if positive:text+='marker\t100\t8 / 8\tsetup_control\n'
                (p/'results_tab.tsv').write_text(text)
        def fake_output(argv,**kwargs):
            argv=list(map(str,argv))
            if '--inspect' in argv:return json.dumps(dict(ready=True,missing=[],environment=str(root/'env'),python=str(root/'env/bin/python'),tools={n:str(root/'env/bin'/n) for n in setup_worker.TOOLS},plasmidfinder_db=str(root/'databases/plasmidfinder'),mob_db=str(root/'databases/mob')))
            if 'list' in argv:return '[]'
            if 'rev-parse' in argv:return 'abc123\n'
            raise AssertionError(argv)
        from contextlib import ExitStack
        with ExitStack() as stack:
            stack.enter_context(patch('setup_worker.bootstrap_conda',return_value=Path('/fake/conda')))
            stack.enter_context(patch('setup_worker.run',side_effect=fake_run))
            stack.enter_context(patch('setup_worker.subprocess.check_output',side_effect=fake_output))
            if hasattr(os,'sched_setaffinity'):stack.enter_context(patch('setup_worker.os.sched_setaffinity'))
            result=setup_worker.install(root)
        return result,calls
    def test_isolated_install_verifies_before_ready(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/'runtime';result,calls=self.exercise(root,True)
            self.assertTrue(result['ready']);self.assertEqual(result['positive_control'],'passed')
            self.assertTrue((root/'ready.json').is_file())
            command,env=calls[0];self.assertIn('--prefix',command);self.assertEqual(command[command.index('--prefix')+1],str(root/'env'))
            self.assertTrue(env['CONDA_PKGS_DIRS'].startswith(str(root)));self.assertTrue(env['TMPDIR'].startswith(str(root)))
            self.assertFalse(any('sudo' in c for c,_ in calls))
            pip_commands=[c for c,_ in calls if '-m' in c and 'pip' in c and 'install' in c]
            self.assertIn('--constraint',pip_commands[0]);self.assertIn('cgecore==1.5.6',pip_commands[0])
            self.assertIn('numpy<1.23.5',(root/'pip-constraints.txt').read_text())
    def test_control_failure_never_writes_ready(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/'runtime'
            with self.assertRaisesRegex(ValueError,'positive control failed'):self.exercise(root,False)
            self.assertFalse((root/'ready.json').exists());self.assertFalse((root/'checks/positive_control.json').exists())
