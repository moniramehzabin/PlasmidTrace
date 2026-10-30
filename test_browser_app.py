import http.client,json,tempfile,threading,time,unittest,zipfile,io
from pathlib import Path
from unittest.mock import patch
from browser_app import State,clean_settings,run_settings,serve

class BrowserAppTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.server=serve(root=self.root);self.thread=threading.Thread(target=self.server.serve_forever,daemon=True);self.thread.start()
        self.prefix='/'+self.server.token+'/'
    def tearDown(self):self.server.shutdown();self.server.server_close();self.thread.join();self.tmp.cleanup()
    def request(self,path='',data=None,headers=None,raw_path=None):
        c=http.client.HTTPConnection('127.0.0.1',self.server.server_port,timeout=20)
        h=dict(headers or {})
        body=None if data is None else json.dumps(data)
        if body:h['Content-Type']='application/json'
        c.request('GET' if data is None else 'POST',raw_path if raw_path is not None else self.prefix+path,body=body,headers=h)
        r=c.getresponse();result=(r.status,dict(r.getheaders()),r.read());c.close();return result
    def test_home_assets_and_private_headers(self):
        status,headers,body=self.request();self.assertEqual(status,200);self.assertIn(b'Find plasmid candidates',body);self.assertIn("frame-ancestors 'none'",headers['Content-Security-Policy'])
        for name in ('app.css','app.js'):
            status,_,body=self.request(name);self.assertEqual(status,200);self.assertGreater(len(body),1000)
    def test_private_token_and_host_required(self):
        self.assertEqual(self.request(raw_path='/')[0],403)
        self.assertEqual(self.request(headers={'Host':'evil.example'})[0],403)
        self.assertEqual(self.request(raw_path=self.prefix+'../CITATION.cff')[0],404)
    def test_cross_origin_post_rejected(self):
        status,_,body=self.request('api/run',{},headers={'Origin':'https://evil.example'});self.assertEqual(status,403)
        self.assertNotIn(b'Traceback',body)
    def test_schema_validation_and_no_arbitrary_settings(self):
        for data in ({'unknown':'x'},{'run_mob':'false'},{'threads':{}},{'backend':'remote'}):
            with self.assertRaises(ValueError):clean_settings(data)
        self.assertEqual(self.request('api/run',{'unknown':'x'})[0],400)
    def test_inactive_settings_do_not_conflict(self):
        s=clean_settings({'workflow':'screen','amrfinder_mode':'run'});self.assertEqual(run_settings(s)['amrfinder_mode'],'off')
        s=clean_settings({'workflow':'cargo','run_barrnap':True,'rrna_gff':'old.gff'});r=run_settings(s);self.assertFalse(r['run_barrnap']);self.assertEqual(r['rrna_gff'],'')
    def test_actual_analysis_and_downloads(self):
        assembly=self.root/'assembly.fa';assembly.write_text('>NODE_1_length_1000_cov_10\n'+'A'*1000+'\n')
        status,_,body=self.request('api/run',{'assembly':str(assembly),'output_root':str(self.root/'results'),'automatic':False,'exploratory':True});self.assertEqual(status,202,body)
        key=json.loads(body)['id'];job=self.server.state.jobs[key];self.assertTrue(job.finished.wait(15))
        status,_,body=self.request('api/status/'+key);d=json.loads(body);self.assertEqual(d['state'],'finished',d['log']);self.assertTrue(d['excel']);self.assertTrue(d['report'])
        status,headers,body=self.request('report/'+key);self.assertEqual(status,200);self.assertIn('sandbox',headers['Content-Security-Policy'])
        status,headers,body=self.request('excel/'+key);self.assertEqual(status,200);self.assertIn('attachment',headers['Content-Disposition'])
        status,_,body=self.request('download/'+key);self.assertEqual(status,200)
        with zipfile.ZipFile(io.BytesIO(body)) as z:self.assertIn('output/evidence.tsv',z.namelist())
        self.assertEqual(json.loads(self.server.state.settings_file.read_text())['assembly'],str(assembly))
    def test_picker_uses_native_path_without_upload(self):
        with patch.object(State,'pick',return_value={'path':r'E:\data\assembly.fa'}):
            status,_,body=self.request('api/pick',{'kind':'file'});self.assertEqual(status,200);self.assertEqual(json.loads(body)['path'],r'E:\data\assembly.fa')
    def test_duplicate_active_run_rejected(self):
        assembly=self.root/'assembly.fa';assembly.write_text('>c\nAAAA\n')
        with patch('browser_app.Job.run'):
            s={'assembly':str(assembly),'output_root':str(self.root/'results'),'automatic':False,'exploratory':True}
            key=self.server.state.start(s)
            self.assertEqual(self.request('api/run',s)[0],400)
            self.server.state.jobs[key].finished.set()
    def test_cancel_endpoint_and_unknown_job(self):
        self.assertEqual(self.request('api/cancel/missing',{})[0],404)
        self.assertEqual(self.request('api/status/missing')[0],404)
    def test_check_reports_missing_selected_tools(self):
        with patch('browser_app.Backend.probe',return_value={'python':'python','tools':{'barrnap':None}}):
            status,_,body=self.request('api/check',{'run_barrnap':True,'automatic':False,'exploratory':True});self.assertEqual(status,200);d=json.loads(body);self.assertFalse(d['ready']);self.assertEqual(d['missing'],['barrnap'])

if __name__=='__main__':unittest.main()
