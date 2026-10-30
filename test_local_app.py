import json, os, subprocess, sys, tempfile, threading, time, unittest, zipfile
from pathlib import Path
from unittest.mock import patch
from local_app import Backend, Job, build_command, cargo_config, make_archive
from app_export import export_results

class LocalAppTests(unittest.TestCase):
    def settings(self, assembly):
        return dict(workflow='screen',assembly=str(assembly),threads='2',memory='7')
    def test_windows_paths_and_shell_quoting(self):
        b=Backend('wsl',environment="env; echo nope")
        self.assertEqual(b.path(r'E:\my folder\reads.fastq.gz'),'/mnt/e/my folder/reads.fastq.gz')
        with tempfile.TemporaryDirectory() as d:
            cmd=build_command(self.settings(r'E:\data\assembly.fa'),b,Path(d))
            self.assertIn("conda activate 'env; echo nope'",cmd[-1])
            self.assertIn('/mnt/e/data/assembly.fa',cmd[-1])
            self.assertEqual(cmd[:3],['wsl.exe','-d','Ubuntu-24.04'])
    def test_windows_wsl_uses_script_file_without_inline_shell_expansion(self):
        import local_app
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            backend = Backend('wsl')
            script = 'pid=$!\nwait "$pid"\nprintf "%s" "$HOME"\n'
            with patch.object(local_app, 'ROOT', root), patch.object(local_app.os, 'name', 'nt'):
                command = backend.shell(script)
            self.assertEqual(command[3:5], ['--exec', 'bash'])
            self.assertNotIn('-lc', command)
            self.assertEqual(Path(command[-1]).read_text(), script)

    def test_conflicting_and_missing_selections(self):
        with tempfile.TemporaryDirectory() as d:
            s=self.settings('a.fa');s.update(run_barrnap=True,rrna_gff='r.gff')
            with self.assertRaisesRegex(ValueError,'barrnap OR'):build_command(s,Backend(),Path(d))
            s=self.settings('a.fa');s.update(workflow='reconstruct')
            with self.assertRaisesRegex(ValueError,'both paired'):build_command(s,Backend(),Path(d))
            s.update(reads1='a',reads2='a')
            with self.assertRaisesRegex(ValueError,'different'):build_command(s,Backend(),Path(d))
    def test_reuse_and_recovery_arguments(self):
        with tempfile.TemporaryDirectory() as d:
            s=self.settings('a.fa');s.update(workflow='reuse',previous='old',reads1='r1',reads2='r2')
            cmd=build_command(s,Backend(),Path(d))
            self.assertIn('--reuse-unicycler',cmd);self.assertNotIn('--recover-from',cmd)
            s['workflow']='recover';cmd=build_command(s,Backend(),Path(d))
            self.assertIn('--recover-from',cmd);self.assertIn('--unicycler-python',cmd)
    def test_cargo_configuration_requires_matching_annotation(self):
        s=dict(geneminer_mode='import',geneminer_result='hits.tsv')
        with self.assertRaisesRegex(ValueError,'genbank'):cargo_config(s,Backend())
        s.update(geneminer_genbank='candidate.gbk',geneminer_workbook='report.xlsx')
        c=cargo_config(s,Backend())['geneminer'];self.assertEqual(c['genbank'],'candidate.gbk');self.assertNotIn('command',c)
    def test_automatic_cargo_commands_are_argument_lists(self):
        s=dict(amrfinder_mode='run',amrfinder_exe='amrfinder',threads='2',virulencefinder_mode='run',virulencefinder_db='/db')
        c=cargo_config(s,Backend());self.assertEqual(c['amrfinder']['command'][1:3],['-n','{assembly}'])
        self.assertIn('{out}',c['virulencefinder']['command'])
    def test_real_screen_job_and_exports(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'assembly.fasta';p.write_text('>NODE_1_length_1000_cov_10\n'+'A'*1000+'\n>NODE_2_length_1000_cov_50\n'+'C'*1000+'\n')
            job=Job(self.settings(p),Backend(),Path(d)/'runs');job.run()
            self.assertEqual(job.state,'finished',(job.folder/'console.log').read_text())
            self.assertTrue((job.folder/'output/evidence.tsv').exists())
            self.assertTrue((job.folder/'report.html').exists())
            self.assertEqual(json.loads((job.folder/'app_run.json').read_text())['state'],'finished')
            with zipfile.ZipFile(job.folder/'PlasmidTrace_results.xlsx') as z:
                self.assertIn('xl/workbook.xml',z.namelist())
            archive=make_archive(job.folder)
            with zipfile.ZipFile(archive) as z:self.assertIn('output/evidence.tsv',z.namelist())
            # New runs do not reuse output directories.
            second=Job(self.settings(p),Backend(),Path(d)/'runs');self.assertNotEqual(second.folder,job.folder)
    def test_failure_is_not_success(self):
        with tempfile.TemporaryDirectory() as d:
            job=Job(self.settings(Path(d)/'missing.fa'),Backend(),Path(d)/'runs');job.run()
            self.assertEqual(job.state,'failed');self.assertNotEqual(json.loads((job.folder/'app_run.json').read_text())['exit_code'],0)
    def test_cancel_before_launch(self):
        with tempfile.TemporaryDirectory() as d:
            job=Job(self.settings('a.fa'),Backend(),Path(d));job.cancel();job.run()
            self.assertEqual(job.state,'cancelled');self.assertIsNone(job.process)
    @unittest.skipIf(os.name=='nt','POSIX process-group test')
    def test_cancel_child_process_group(self):
        with tempfile.TemporaryDirectory() as d:
            job=Job(self.settings('a.fa'),Backend(),Path(d));job.command=[sys.executable,'-u','-c','import time;print("started",flush=True);time.sleep(60)']
            t=threading.Thread(target=job.run);t.start()
            for _ in range(100):
                if job.state=='running':break
                time.sleep(.02)
            job.cancel();t.join(10);self.assertFalse(t.is_alive());self.assertEqual(job.state,'cancelled')
    def test_report_escapes_input_and_workbook_has_no_formulas(self):
        from local_app import export_report
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'output').mkdir();(root/'output/report.md').write_text('<script>alert(1)</script>')
            (root/'output/table.tsv').write_text('contig\texplanation\n=HYPERLINK("evil")\t<script>\n')
            export_report(root,'finished');text=(root/'report.html').read_text();self.assertNotIn('<script>',text)
            with zipfile.ZipFile(root/'PlasmidTrace_results.xlsx') as z:
                for n in z.namelist():
                    if n.startswith('xl/worksheets/'):self.assertNotIn('<f>',z.read(n).decode())
    @unittest.skipIf(os.name=='nt','Linux wrapper integration test')
    def test_wsl_shell_wrapper_executes_screening_with_space_paths(self):
        with tempfile.TemporaryDirectory(prefix='pt space ') as d:
            p=Path(d)/'input fasta.fa';p.write_text('>contig\n'+'A'*1000+'\n')
            b=Backend('wsl');b.shell=lambda script:['bash','-lc',script]
            job=Job(self.settings(p),b,Path(d)/'results');job.run()
            self.assertEqual(job.state,'finished',(job.folder/'console.log').read_text())
            self.assertTrue((job.folder/'linux_process_group.pid').is_file())
    @unittest.skipIf(os.name=='nt','Linux wrapper integration test')
    def test_wsl_wrapper_cancel_marker_prevents_launch(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'cancel.requested').touch()
            b=Backend('wsl');b.shell=lambda script:['bash','-lc',script]
            result=subprocess.run(b.wrap(['-c','print("should not run")'],root),capture_output=True,text=True)
            self.assertEqual(result.returncode,130);self.assertNotIn('should not run',result.stdout)
    def test_missing_environment_tools_stops_early(self):
        with tempfile.TemporaryDirectory() as d:
            s=self.settings('a.fa');s['run_barrnap']=True
            with patch.object(Backend,'probe',return_value={'tools':{'barrnap':None}}):
                job=Job(s,Backend(),Path(d));job.run()
            self.assertEqual(job.state,'failed');self.assertIsNone(job.process)
            self.assertIn('Missing tools',(job.folder/'console.log').read_text())

if __name__=='__main__':unittest.main()
