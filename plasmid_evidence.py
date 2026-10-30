#!/usr/bin/env python3
"""Conservative plasmid evidence triage. Python 3.9+, no Python dependencies."""
import argparse, csv, hashlib, json, re, shutil, subprocess, sys
from pathlib import Path
from statistics import median
import finder_adapter as finder
import annotation_evidence as annotation
VERSION = '1.0.0'

def fasta(path):
    records = {}; name = None
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split()[0]
            if name in records: raise ValueError('Duplicate FASTA ID: '+name)
            records[name] = ''
        elif line:
            if name is None: raise ValueError('Sequence before FASTA header')
            records[name] += line.upper()
    if not records or any(not s or set(s)-set('ACGTRYSWKMBDHVN') for s in records.values()):
        raise ValueError('Empty or invalid nucleotide FASTA')
    return records

def rrna_intervals(path, records):
    intervals = {n: [] for n in records}
    for line in Path(path).read_text().splitlines():
        if line.startswith('#') or not line.strip(): continue
        f = line.split('\t')
        if len(f) != 9: raise ValueError('Invalid GFF row')
        if f[2].lower() != 'rrna': continue
        if f[0] not in records: raise ValueError('rRNA GFF ID not in assembly: '+f[0])
        start, end = int(f[3])-1, int(f[4])
        if not 0 <= start < end <= len(records[f[0]]): raise ValueError('Invalid rRNA coordinates')
        intervals[f[0]].append((start,end))
    fractions = {}
    for n, regions in intervals.items():
        covered = 0; last = -1
        for start,end in sorted(regions):
            covered += max(0, end-max(start,last)); last = max(last,end)
        fractions[n] = covered/len(records[n])
    return fractions

def classify(mob, rrna, elevated, circular):
    reasons = []
    if rrna is not None and rrna >= .5:
        return 'rRNA-dominated warning', 'Most sequence overlaps predicted rRNA; do not accept a plasmid call from this alone.'
    if mob == 'plasmid': reasons.append('MOB-Recon predicts plasmid')
    if elevated: reasons.append('Coverage exceeds exploratory baseline threshold')
    if circular: reasons.append('Same-orientation GFA self-link; closure is not read-validated')
    if reasons: return 'candidate for review', '; '.join(reasons)
    return 'no positive evidence in supplied checks', 'Not evidence of plasmid absence; missing checks remain unassessed.'

def blast_import(paths, records, aliases=None, tsv12=False):
    aliases = aliases or {}
    """Read saved NCBI pairwise text or the documented eight-column TSV."""
    hits = {n: [] for n in records}
    for path in paths:
        text = path.read_text()
        if '<html' in text.lower() or '<?xml' in text.lower():
            raise ValueError('Use plain-text pairwise BLAST or the documented TSV, not HTML/XML')
        if 'Query=' in text:
            query = None; accession = None; title = ''; collecting = False; current = None
            for line in text.splitlines():
                if line.startswith('Query='):
                    query = line.split('=',1)[1].strip().split()[0]
                    query = aliases.get(query,query)
                    if query not in records: raise ValueError('BLAST query ID mismatch: '+query)
                    accession = None; current = None
                elif line.startswith('>'):
                    current = None
                    parts = line[1:].strip().split(None,1)
                    accession = parts[0]; title = parts[1] if len(parts)>1 else ''; collecting = True
                elif collecting and line.startswith('Length='):
                    collecting = False
                elif collecting:
                    title += ' '+line.strip()
                elif accession and 'Identities =' in line:
                    m = re.search(r'Identities\s*=\s*(\d+)/(\d+)',line)
                    if not m: raise ValueError('Cannot parse BLAST identities')
                    identity = 100*int(m[1])/int(m[2])
                    current = dict(accession=accession,title=title.strip(),identity=identity,qcov=None,alignment_length=int(m[2]),query_intervals=[])
                    hits[query].append(current)
                elif current is not None:
                    m = re.match(r'^Query\s+(\d+)\s+[A-Za-z-]+\s+(\d+)\s*$',line)
                    if m:
                        lo, hi = sorted((int(m[1]),int(m[2])))
                        if not 1 <= lo <= hi <= len(records[query]): raise ValueError('BLAST query coordinates outside assembly sequence')
                        current['query_intervals'].append([lo,hi])
        else:
            for line in text.splitlines():
                if not line.strip() or line.startswith('#'): continue
                f = line.split('\t')
                f[0] = aliases.get(f[0],f[0])
                if len(f)==12 and tsv12:
                    # Explicit schema from this project's saved results; not default BLAST outfmt 6.
                    identity, length, qlen = float(f[2]), int(f[3]), int(f[4])
                    lo,hi = sorted((int(f[6]),int(f[7])))
                    if f[0] not in records: raise ValueError('BLAST query ID mismatch: '+f[0])
                    if qlen!=len(records[f[0]]) or not 1<=lo<=hi<=qlen or not 0<=identity<=100: raise ValueError('Invalid 12-column BLAST values')
                    hits[f[0]].append(dict(accession=f[1],title='Title unavailable in imported table',identity=identity,qcov=None,alignment_length=length,query_intervals=[[lo,hi]]))
                    continue
                if len(f)!=8: raise ValueError('Expected documented eight- or twelve-column BLAST TSV')
                if f[0] not in records: raise ValueError('BLAST query ID mismatch: '+f[0])
                identity, qcov = float(f[2]), float(f[4])
                if not 0<=identity<=100 or not 0<=qcov<=100: raise ValueError('Invalid BLAST percentage')
                hits[f[0]].append(dict(accession=f[1],title=f[7],identity=identity,qcov=qcov,alignment_length=int(f[3]),query_intervals=[]))
    return hits

def blast_summary(hits):
    if not hits: return 'no imported hits; search completion not verified'
    kinds = set()
    for h in hits:
        title = h['title'].lower()
        kinds.add('plasmid' if 'plasmid' in title else 'chromosome' if 'chromosome' in title else 'unspecified')
    return ', '.join(sorted(kinds))+' record matches; location in this isolate remains unresolved'

def reference_summaries(hits, query_length):
    grouped = {}
    for h in hits: grouped.setdefault(h['accession'],[]).append(h)
    rows = []
    for accession, hs in grouped.items():
        regions = sorted({tuple(r) for h in hs for r in h['query_intervals']})
        covered = 0; last = 0
        for lo,hi in regions:
            covered += max(0,hi-max(lo-1,last)); last = max(last,hi)
        complete_coords = all(h['query_intervals'] for h in hs)
        coverage = 100*covered/query_length if complete_coords else None
        if coverage is None:
            reported = [h['qcov'] for h in hs if h['qcov'] is not None]
            coverage = max(reported) if reported else None
        longest = max(hs,key=lambda h:h['alignment_length'])
        rows.append(dict(accession=accession,title=hs[0]['title'],query_coverage_percent=coverage,
                         coverage_source='coordinate union' if complete_coords else 'reported qcovs' if coverage is not None else 'unavailable',
                         alignment_sections=len(hs),longest_alignment_bp=longest['alignment_length'],
                         longest_alignment_identity_percent=longest['identity']))
    return rows

def exact_aliases(path, records):
    aliases = {}
    complement = str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')
    for name, sequence in fasta(path).items():
        reverse = sequence.translate(complement)[::-1]
        matches = [n for n,s in records.items() if s==sequence or s==reverse]
        if len(matches)!=1: raise ValueError('Alias sequence must match exactly one assembly record: '+name)
        if name in records and name!=matches[0]: raise ValueError('Alias conflicts with assembly ID: '+name)
        aliases[name] = matches[0]
    return aliases

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='version', version=VERSION)
    p.add_argument('--assembly', type=Path)
    p.add_argument('--unicycler-dir', type=Path, help='Reuse assembly.fasta and assembly.gfa; does not run Unicycler')
    p.add_argument('--blast-tsv12', action='store_true', help='Explicit custom schema: qseqid sacc pident length qlen slen qstart qend sstart send evalue bitscore')
    p.add_argument('--alias-fasta', type=Path, help='Map renamed BLAST queries by unique exact DNA or reverse-complement match')
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--rrna-gff', type=Path)
    p.add_argument('--mob-report', type=Path)
    p.add_argument('--blast-results', type=Path, nargs='+', default=[], help='Saved pairwise BLAST text or eight-column TSV; no network calls')
    p.add_argument('--gfa', type=Path, help='Graph IDs must match FASTA IDs; no guessed mapping')
    p.add_argument('--run-barrnap', action='store_true')
    p.add_argument('--run-mob', action='store_true', help='Runs installed MOB-Recon; its database may require download')
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--coverage-ratio', type=float, default=2.0, help='Exploratory cutoff, not validated plasmid score')
    p.add_argument('--annotation-genbank',type=Path,help='PGAP/Prokka nucleotide GenBank, mapped by exact DNA')
    p.add_argument('--annotation-gff',type=Path,help='PGAP/Prokka GFF3 gene evidence; not a plasmid classifier')
    p.add_argument('--annotation-fasta',type=Path,help='Annotation-associated genomic FASTA for exact same-DNA ID mapping')
    finder.arguments(p)
    a = p.parse_args()
    finder.validate(a)
    if a.threads < 1 or a.coverage_ratio <= 0: p.error('Threads and coverage ratio must be positive')
    if a.run_barrnap and a.rrna_gff: p.error('Use --run-barrnap OR --rrna-gff')
    if a.run_mob and a.mob_report: p.error('Use --run-mob OR --mob-report')
    if a.unicycler_dir:
        if a.assembly or a.gfa: p.error('Use --unicycler-dir OR --assembly/--gfa')
        a.assembly = a.unicycler_dir/'assembly.fasta'
        a.gfa = a.unicycler_dir/'assembly.gfa'
    if not a.assembly: p.error('Provide --assembly or --unicycler-dir')
    records = fasta(a.assembly)
    if a.annotation_fasta and not a.annotation_gff: p.error('--annotation-fasta requires --annotation-gff')
    if a.annotation_genbank and a.annotation_gff: p.error('Choose --annotation-genbank OR --annotation-gff')
    annotation_aliases = exact_aliases(a.annotation_fasta,records) if a.annotation_fasta else {}
    genes,declared = annotation.read_gff(a.annotation_gff,records,annotation_aliases) if a.annotation_gff else ({n:[] for n in records},{n:[] for n in records})
    if a.annotation_genbank: genes,declared,annotation_aliases = annotation.read_genbank(a.annotation_genbank,records)
    annotation_checked = bool(a.annotation_gff or a.annotation_genbank)
    aliases = exact_aliases(a.alias_fasta,records) if a.alias_fasta else {}
    blast = blast_import(a.blast_results, records,aliases,a.blast_tsv12)
    if a.out.exists(): raise ValueError('Output folder exists; choose a new folder to preserve results')
    a.out.mkdir(parents=True)
    (a.out/'annotation_evidence.json').write_text(json.dumps({'genes':genes,'source_plasmid_labels':declared,'aliases':annotation_aliases},indent=2)+'\n')
    commands = []
    def run(cmd, log, stdout=None):
        commands.append(cmd); print('Running:', ' '.join(cmd), flush=True)
        if not shutil.which(cmd[0]): raise ValueError('Missing executable: '+cmd[0])
        with log.open('w') as err:
            if stdout:
                with stdout.open('w') as out: subprocess.run(cmd,stdout=out,stderr=err,check=True)
            else: subprocess.run(cmd,stdout=err,stderr=err,check=True)
    if a.run_barrnap:
        a.rrna_gff = a.out/'rrna.gff'
        run(['barrnap','--kingdom','bac','--threads',str(a.threads),str(a.assembly.resolve())],a.out/'barrnap.log',a.rrna_gff)
    if a.run_mob:
        run([a.mob_executable,'-i',str(a.assembly.resolve()),'-o',str((a.out/'mob_recon').resolve()),'-n',str(a.threads)],a.out/'mob_recon.log')
        a.mob_report = a.out/'mob_recon'/'contig_report.txt'
    if a.run_plasmidfinder:
        run(finder.command(a,a.assembly,a.out/'plasmidfinder'),a.out/'plasmidfinder.log')
        a.plasmidfinder_results = a.out/'plasmidfinder/results_tab.tsv'
    markers = finder.import_hits(a.plasmidfinder_results,records,a.plasmidfinder_identity,a.plasmidfinder_coverage) if a.plasmidfinder_results else {n:[] for n in records}
    (a.out/'plasmidfinder_hits.json').write_text(json.dumps(markers,indent=2)+'\n')
    rrna = rrna_intervals(a.rrna_gff,records) if a.rrna_gff else {}
    mob = {}
    if a.mob_report:
        with a.mob_report.open() as handle:
            reader = csv.DictReader(handle,delimiter='\t')
            if not {'contig_id','molecule_type'} <= set(reader.fieldnames or []): raise ValueError('MOB report lacks required columns')
            for row in reader:
                n = row['contig_id']
                if n not in records: raise ValueError('MOB ID mismatch: '+n)
                mob[n] = row['molecule_type']
    links = set(); graph_ids = set(); external_neighbors = {}
    if a.gfa:
        for line in a.gfa.read_text().splitlines():
            f = line.split('\t')
            if f[0] == 'S': graph_ids.add(f[1])
            if f[0] == 'L' and len(f)>=6:
                if f[1]==f[3] and f[2]==f[4]: links.add(f[1])
                if f[1]!=f[3]:
                    external_neighbors.setdefault(f[1],set()).add(f[3])
                    external_neighbors.setdefault(f[3],set()).add(f[1])
    # Validate same-ID graph DNA before assigning a self-link to an assembly record.
    if a.gfa:
        for line in a.gfa.read_text().splitlines():
            f = line.split('\t')
            if f[0]=='S' and f[1] in records and f[2]!='*' and f[2].upper()!=records[f[1]]:
                raise ValueError('GFA sequence differs from assembly: '+f[1])
    coverage = {}
    for n in records:
        m = re.search(r'_cov_([0-9]+(?:\.[0-9]+)?)',n)
        if m: coverage[n] = float(m[1])
    baseline_values = [coverage[n] for n in records if n in coverage and len(records[n])>=10000]
    baseline = median(baseline_values) if baseline_values else None
    fields = ['contig','length_bp','gc_percent','coverage_header','coverage_ratio','rrna_fraction','mob_classification','gfa_self_link','status','explanation','blast_evidence','graph_context','external_graph_neighbors']
    fields += ['plasmidfinder_markers','annotation_gene_evidence','annotation_plasmid_labels']
    rows = []; candidates = []
    for n,s in records.items():
        ratio = coverage[n]/baseline if baseline and n in coverage else None
        connected_loop = n in links and bool(external_neighbors.get(n))
        circular = n in links and not connected_loop
        status, explanation = classify(mob.get(n),rrna.get(n),ratio is not None and ratio>=a.coverage_ratio,circular)
        if genes[n]:
            explanation += '; annotation evidence: '+', '.join(h['product'] or h['gene'] for h in genes[n])
            if status != 'rRNA-dominated warning' and annotation.nominates(genes[n]): status = 'candidate for review'
            elif status == 'no positive evidence in supplied checks': status = 'annotation context for review'
        if declared[n]:
            explanation += '; source annotation plasmid label (supplied metadata, not independent detection): '+', '.join(declared[n])
        if markers[n]:
            explanation += '; PlasmidFinder marker(s): '+', '.join(h['marker'] for h in markers[n])
            if status != 'rRNA-dominated warning': status = 'candidate for review'
        if connected_loop:
            if status == 'no positive evidence in supplied checks':
                status = 'connected graph loop warning'
            explanation += '; self-link also connects to other segments; not an isolated circular sequence'
        elif circular and status == 'candidate for review':
            explanation += '; isolated single-segment circular candidate, not a confirmed plasmid'
        if blast[n]:
            explanation += '; BLAST: '+blast_summary(blast[n])
        if status == 'candidate for review': candidates.append(n)
        gc = 100*(s.count('G')+s.count('C'))/sum(s.count(b) for b in 'ACGT') if any(b in s for b in 'ACGT') else None
        rows.append(dict(zip(fields,[n,len(s),gc,coverage.get(n,'not assessed'),ratio if ratio is not None else 'not assessed',rrna.get(n,'not assessed'),mob.get(n,'not assessed'),'yes' if n in links else ('no self-link detected' if n in graph_ids else 'not assessed'),status,explanation,blast_summary(blast[n]) if a.blast_results else 'not assessed','connected self-loop' if connected_loop else 'isolated same-orientation self-loop' if circular else 'no same-orientation self-loop detected' if n in graph_ids else 'not assessed',','.join(sorted(external_neighbors.get(n,set())))])))
    for row in rows:
        row['annotation_gene_evidence'] = json.dumps(genes[row['contig']]) if annotation_checked else 'not assessed'
        row['annotation_plasmid_labels'] = ','.join(declared[row['contig']]) if annotation_checked else 'not assessed'
        row['plasmidfinder_markers'] = ','.join(h['marker'] for h in markers[row['contig']]) if a.plasmidfinder_results else 'not assessed'
    with (a.out/'evidence.tsv').open('w') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    with (a.out/'candidates.fasta').open('w') as h:
        for n in candidates: h.write('>'+n+'\n'+records[n]+'\n')
    inputs = {k:str(v.resolve()) for k,v in [('assembly',a.assembly),('rrna_gff',a.rrna_gff),('mob_report',a.mob_report),('gfa',a.gfa)] if v}
    if a.annotation_genbank: inputs['annotation_genbank'] = str(a.annotation_genbank.resolve())
    if a.annotation_gff: inputs['annotation_gff'] = str(a.annotation_gff.resolve())
    if a.annotation_fasta: inputs['annotation_fasta'] = str(a.annotation_fasta.resolve())
    if a.plasmidfinder_results: inputs['plasmidfinder_results'] = str(a.plasmidfinder_results.resolve())
    if a.alias_fasta: inputs['alias_fasta'] = str(a.alias_fasta.resolve())
    (a.out/'sequence_aliases.json').write_text(json.dumps(aliases,indent=2)+'\n')
    inputs.update({f'blast_{i+1}':str(path.resolve()) for i,path in enumerate(a.blast_results)})
    (a.out/'blast_hits.json').write_text(json.dumps(blast,indent=2)+'\n')
    checksums = {k:hashlib.sha256(Path(v).read_bytes()).hexdigest() for k,v in inputs.items()}
    (a.out/'provenance.json').write_text(json.dumps({'version':VERSION,'finder_database_sha256':finder.database_fingerprint(a.plasmidfinder_db),'plasmidfinder_identity':a.plasmidfinder_identity,'plasmidfinder_coverage':a.plasmidfinder_coverage,'inputs':inputs,'sha256':checksums,'commands':commands,'coverage_baseline':baseline,'baseline_definition':'median SPAdes header coverage of contigs >=10000 bp','coverage_ratio_cutoff':a.coverage_ratio},indent=2)+'\n')
    warnings = ['No result is a confirmed plasmid. This is evidence triage, not a validated classifier.', 'Coverage is original SPAdes header coverage, not mapped base depth or copy number.', 'Graph support distinguishes isolated and externally connected matching-ID same-orientation self-links; multi-segment circles and SPAdes numeric-ID mapping are not assessed.', 'BLAST imports retain record titles and identity; chromosome matches do not exclude plasmids. Pairwise query coverage is calculated from the union of aligned query coordinates; this includes all imported matching sections regardless of identity. Missing coordinates remain unassessed. Identity shown is for the longest matching section, not a whole-reference average. Gene interpretation, reassembly and junction read validation are not implemented.', 'Absence of candidates does not establish plasmid absence.']
    if a.gfa and not (graph_ids & records.keys()): warnings.append('WARNING: GFA and FASTA IDs do not match. Graph checks are unassessed for every contig.')
    report = ['# PlasmidTrace v'+VERSION, f'Contigs: {len(records)}; review candidates: {len(candidates)}', '## Limits']+['- '+w for w in warnings]
    summaries = {n:reference_summaries(hs,len(records[n])) for n,hs in blast.items() if hs}
    (a.out/'blast_reference_summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    if a.unicycler_dir: report.append('Existing Unicycler assembly and graph reused; no assembly was run. Circular self-links remain unvalidated by reads.')
    if aliases: report.append('Renamed query IDs mapped by unique exact DNA or reverse-complement matches: '+str(aliases))
    report += ['## Finder checks', 'MOB-Recon: '+('report assessed' if a.mob_report else 'not run or supplied'), 'PlasmidFinder: '+(str(sum(map(len,markers.values())))+' qualifying marker hits; coverage refers to marker length' if a.plasmidfinder_results else 'not run or supplied'), 'Marker absence does not exclude an unrepresented plasmid. rRNA warnings take priority and retain marker evidence for review.']
    report += ['Annotation check: '+('Annotation product/gene keyword evidence assessed; no protein homology search was run' if annotation_checked else 'not assessed'), 'Replication/partition annotations nominate review candidates, not confirmed plasmids. Generic partition and mobility genes are contextual evidence because they may occur on chromosomes. Source plasmid labels are supplied metadata. Noncircular candidates remain in evidence.tsv even if junction testing is skipped.']
    report += ['## Results']+[f'- {r["contig"]}: {r["status"]}. {r["explanation"]}' for r in rows if r['status']!='no positive evidence in supplied checks']
    for n, refs in summaries.items():
        report.append('### BLAST reference matches: '+n)
        for r in refs:
            cov = f"{r['query_coverage_percent']:.2f}%" if r['query_coverage_percent'] is not None else 'unavailable'
            report.append(f"- {r['accession']}: query coverage {cov} ({r['coverage_source']}); longest alignment {r['longest_alignment_bp']} bp at {r['longest_alignment_identity_percent']:.3f}% identity; {r['alignment_sections']} matching sections. {r['title']}")
    (a.out/'report.md').write_text('\n\n'.join(report)+'\n')
    print('Completed. Read',a.out/'report.md')
if __name__ == '__main__':
    try: main()
    except (ValueError,OSError,subprocess.CalledProcessError) as e:
        print('ERROR:',e,file=sys.stderr);sys.exit(1)
