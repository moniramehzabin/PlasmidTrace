import json, tempfile, unittest
from pathlib import Path
from cargo_analysis import normalize, analyze

class CargoTests(unittest.TestCase):
    def setUp(self):
        self.t=tempfile.TemporaryDirectory();self.root=Path(self.t.name)
    def tearDown(self): self.t.cleanup()
    def table(self,text):
        p=self.root/'data.tsv';p.write_text(text);return p
    def test_amr_preserves_element_type(self):
        p=self.table('Contig id\tGene symbol\tElement type\tMethod\tStart\tStop\nc\tstressX\tSTRESS\tBLAST\t2\t8\n')
        self.assertEqual(normalize('amrfinder',p,{'c':'A'*10},{})[0]['category'],'STRESS')
    def test_amr_current_schema(self):
        p=self.table('Contig id\tElement symbol\tType\tMethod\tStart\tStop\nc\tg\tAMR\tBLAST\t2\t8\n')
        self.assertEqual(normalize('amrfinder',p,{'c':'A'*10},{})[0]['gene'],'g')
    def test_unknown_contig_fails(self):
        p=self.table('Contig id\tGene symbol\tElement type\tMethod\tStart\tStop\nx\tg\tAMR\tBLAST\t2\t8\n')
        with self.assertRaises(ValueError): normalize('amrfinder',p,{'c':'A'*10},{})
    def test_geneminer_rejected_row_preserved(self):
        p=self.table('Locus_tag\tFamily_target\tEvidence_class\tConfidence\tDecision\nL1\tmerA\tSupporting\tWeak\tReject\n')
        hit=normalize('geneminer',p,{'c':'A'*10},{'L1':'c'})[0]
        self.assertEqual(hit['decision'],'Reject')
        with self.assertRaises(ValueError): normalize('geneminer',p,{'c':'A'*10},{})
    def test_virulence_reverse_coordinates(self):
        p=self.table('Virulence factor\tIdentity\tContig\tPosition in contig\nv\t99\tc\t8..2\n')
        self.assertEqual(normalize('virulencefinder',p,{'c':'A'*10},{})[0]['start'],2)
    def test_import_empty_is_assessed_zero(self):
        fa=self.root/'a.fa';fa.write_text('>c\nAAAAAAAAAA\n')
        p=self.table('Contig id\tElement symbol\tType\tMethod\tStart\tStop\n')
        config=self.root/'cfg.json';config.write_text(json.dumps({'amrfinder':{'result':str(p)}}))
        result=analyze(fa,config,self.root/'out')
        self.assertEqual(result['modules']['amrfinder']['matched_rows'],0)
        self.assertEqual((self.root/'out/amrfinder/imported_results.tsv').read_bytes(),p.read_bytes())
        with self.assertRaises(ValueError): analyze(fa,config,self.root/'out')
    def test_empty_file_is_not_negative(self):
        with self.assertRaises(ValueError): normalize('amrfinder',self.table(''),{'c':'A'*10},{})

if __name__=='__main__': unittest.main()

class EmptyWorkbookTests(unittest.TestCase):
    def test_empty_and_populated(self):
        import zipfile
        from cargo_analysis import verify_empty_workbook
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'report.xlsx'
            for populated in (False,True):
                with zipfile.ZipFile(path,'w') as z:
                    z.writestr('xl/workbook.xml','<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="All_candidates" r:id="r1"/></sheets></workbook>')
                    z.writestr('xl/_rels/workbook.xml.rels','<Relationships><Relationship Id="r1" Target="worksheets/sheet1.xml"/></Relationships>')
                    cells='<row r="1"><c r="A1" t="inlineStr"><is><t>Gene</t></is></c></row>' if populated else ''
                    z.writestr('xl/worksheets/sheet1.xml','<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>'+cells+'</sheetData></worksheet>')
                if populated:
                    with self.assertRaises(ValueError): verify_empty_workbook(path)
                else: self.assertIn('verified',verify_empty_workbook(path)['status'])
