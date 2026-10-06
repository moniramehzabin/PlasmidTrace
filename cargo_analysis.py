#!/usr/bin/env python3
"""Optional functional cargo adapters. Cargo never changes plasmid classification."""
import argparse, csv, hashlib, json, re, subprocess, zipfile, posixpath
import xml.etree.ElementTree as ET
from pathlib import Path
from plasmid_evidence import fasta
from annotation_evidence import read_genbank
VERSION = '1.0.0'

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''): h.update(block)
    return h.hexdigest()

def verify_empty_workbook(path):
    """Require a real, empty All_candidates sheet; never infer zero from blank TSV."""
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    relns = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'
    with zipfile.ZipFile(path) as z:
        book = ET.fromstring(z.read('xl/workbook.xml'))
        sheets = [x for x in book.findall('s:sheets/s:sheet',ns) if x.get('name') == 'All_candidates']
        if len(sheets) != 1: raise ValueError('Workbook must contain one All_candidates sheet')
        rels = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        targets = [x.get('Target') for x in rels if x.get('Id') == sheets[0].get(relns) and x.get('TargetMode') != 'External']
        if len(targets) != 1: raise ValueError('Invalid workbook sheet relationship')
        target = targets[0]
        member = target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target)
        sheet = ET.fromstring(z.read(member))
        for cell in sheet.findall('.//s:sheetData/s:row/s:c',ns):
            if cell.find('s:v',ns) is not None or cell.find('s:f',ns) is not None or cell.find('s:is',ns) is not None:
                raise ValueError('Workbook All_candidates is populated; blank TSV is inconsistent')
    return {'status':'zero candidates verified in empty All_candidates sheet','workbook_sha256':digest(path)}

def rows(path, required):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames or not set(required) <= set(reader.fieldnames):
            raise ValueError('Unsupported or empty table: '+str(path)+'; expected '+', '.join(required))
        if set(required) == {'Contig id','Method','Start','Stop'} and not (
            {'Element symbol','Type'} <= set(reader.fieldnames) or
            {'Gene symbol','Element type'} <= set(reader.fieldnames)):
            raise ValueError('Unsupported AMRFinder element columns')
        result = list(reader)
        if any(None in r or any(v is None for v in r.values()) for r in result):
            raise ValueError('Malformed TSV row: '+str(path))
        return result

def locus_map(genbank, records):
    # Verify every annotation record against actual assembly DNA before attaching loci.
    _, _, aliases = read_genbank(Path(genbank), records)
    result = {}
    for block in Path(genbank).read_text().split('\n//'):
        if not block.strip(): continue
        name = aliases[re.search(r'^LOCUS\s+(\S+)', block, re.M)[1]]
        for feature in re.split(r'(?m)^     (?=\S)', block.split('FEATURES',1)[-1].split('\nORIGIN',1)[0])[1:]:
            if not feature.startswith('CDS '): continue
            m = re.search(r'/locus_tag="([^"]+)"', feature)
            if m:
                if m[1] in result: raise ValueError('Duplicate annotation locus: '+m[1])
                result[m[1]] = name
    return result

def normalize(tool, raw, records, loci):
    if tool == 'amrfinder':
        required = ['Contig id','Method','Start','Stop']
    elif tool == 'virulencefinder':
        required = ['Virulence factor','Identity','Contig','Position in contig']
    else:
        required = ['Locus_tag','Family_target','Evidence_class','Confidence','Decision']
    output = []
    for r in rows(raw, required):
        if tool == 'amrfinder':
            if 'Element symbol' in r and 'Type' in r:
                gene, category = r['Element symbol'], r['Type']
            elif 'Gene symbol' in r and 'Element type' in r:
                gene, category = r['Gene symbol'], r['Element type']
            else:
                raise ValueError('Unsupported AMRFinder element columns')
            contig, confidence = r['Contig id'], r['Method']
            start, end = sorted((int(r['Start']),int(r['Stop'])))
        elif tool == 'virulencefinder':
            contig, gene = r['Contig'],r['Virulence factor']
            category, confidence = 'virulence-associated','database match'
            coords = re.fullmatch(r'(\d+)\.\.(\d+)',r['Position in contig'])
            if not coords: raise ValueError('Unsupported VirulenceFinder coordinates')
            start,end = sorted(map(int,coords.groups()))
            identity = float(r['Identity'])
            if not 0 <= identity <= 100: raise ValueError('Invalid identity')
        else:
            contig = loci.get(r['Locus_tag'],'')
            gene,category,confidence = r['Family_target'],r['Evidence_class'],r['Confidence']
            start=end=None
            if not contig: raise ValueError('Gene Miner locus not found in verified GenBank: '+r['Locus_tag'])
        if contig not in records: raise ValueError('Cargo contig missing from assembly: '+contig)
        if start is not None and not 1 <= start <= end <= len(records[contig]):
            raise ValueError('Cargo coordinates outside assembly: '+contig)
        output.append(dict(tool=tool,contig=contig,gene=gene,category=category,
            confidence_or_method=confidence,start=start,end=end,
            decision=r.get('Decision','reported match'),raw=r))
    return output

def analyze(assembly, config, out):
    records = fasta(Path(assembly))
    out = Path(out)
    if out.exists() and any(out.iterdir()): raise ValueError('Cargo output must be new or empty')
    out.mkdir(parents=True,exist_ok=True)
    cfg = json.loads(Path(config).read_text())
    if not isinstance(cfg,dict) or set(cfg)-{'amrfinder','virulencefinder','geneminer'} or not cfg:
        raise ValueError('Config requires one or more named cargo modules')
    hits=[]; manifest=dict(version=VERSION,assembly_sha256=digest(assembly),config_sha256=digest(config),modules={})
    for tool, settings in cfg.items():
        if not isinstance(settings,dict) or set(settings)-{'command','result','genbank','provenance_files','workbook'}:
            raise ValueError('Unsupported module settings: '+tool)
        rawdir=out/tool;rawdir.mkdir()
        substitutions={'assembly':str(Path(assembly).resolve()),'out':str(rawdir.resolve())}
        command = settings.get('command')
        if command:
            if not isinstance(command,list) or not all(isinstance(x,str) for x in command):
                raise ValueError('Command must be an argv list; shell commands are not supported')
            argv=[x.replace('{assembly}',substitutions['assembly']).replace('{out}',substitutions['out']) for x in command]
            print('Running cargo module:',tool,flush=True)
            with (rawdir/'command.log').open('wb') as log:
                subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,check=True)
        else: argv=None
        path=Path(settings['result'].replace('{out}',substitutions['out']))
        loci={}
        if tool=='geneminer':
            if not settings.get('genbank'): raise ValueError('Gene Miner requires its original GenBank annotation')
            loci=locus_map(settings['genbank'],records)
        empty_verification=None
        if tool == 'geneminer' and not path.read_text().strip():
            workbook_setting=settings.get('workbook')
            if not workbook_setting: raise ValueError('Blank Gene Miner TSV requires explicit workbook setting')
            workbook=Path(workbook_setting.replace('{out}',substitutions['out']))
            empty_verification=verify_empty_workbook(workbook)
            normalized=[]
        else:
            normalized=normalize(tool,path,records,loci)
        # Preserve original output bytes even in import mode.
        (rawdir/'imported_results.tsv').write_bytes(path.read_bytes())
        provenance={str(Path(p).resolve()):digest(p) for p in settings.get('provenance_files',[])}
        if empty_verification: provenance[str(workbook.resolve())]=digest(workbook)
        if settings.get('genbank'): provenance[str(Path(settings['genbank']).resolve())]=digest(settings['genbank'])
        manifest['modules'][tool]=dict(command=argv,result=str(path.resolve()),result_sha256=digest(path),
            matched_rows=len(normalized),provenance_files=provenance,empty_result_verification=empty_verification,
            status='completed; no reported matches' if not normalized else 'completed; reported matches')
        hits.extend(normalized)
    manifest['state']='completed'
    (out/'cargo_hits.json').write_text(json.dumps(hits,indent=2)+'\n')
    (out/'provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
    columns=['tool','contig','gene','category','confidence_or_method','start','end','decision']
    with (out/'cargo.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(hits)
    text=['# Functional cargo analysis','Cargo results do not change plasmid status. Contig location is an assembly assignment, not proof of plasmid carriage.']
    text += [f'{tool}: {m["matched_rows"]} reported rows.' for tool,m in manifest['modules'].items()]
    text += ['AMRFinderPlus element types are retained separately; stress and virulence entries are not automatically AMR genes.',
        'VirulenceFinder database coverage depends on species. No hit does not establish absence of virulence.',
        'Gene Miner review/reject decisions and confidence are preserved. Predictions do not establish resistance, pathogenicity or bioremediation activity.']
    (out/'report.md').write_text('\n\n'.join(text)+'\n')
    return manifest

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version',action='version',version=VERSION)
    p.add_argument('--assembly',required=True,type=Path)
    p.add_argument('--config',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    a=p.parse_args()
    try: analyze(a.assembly,a.config,a.out)
    except (ValueError,KeyError,OSError,subprocess.CalledProcessError) as e: p.exit(1,'ERROR: '+str(e)+'\n')
    print('Completed. Read',a.out/'report.md')
if __name__=='__main__': main()
