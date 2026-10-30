#!/usr/bin/env python3
"""Active short-read plasmid reconstruction and junction support workflow."""
import argparse, csv, hashlib, json, os, re, shutil, subprocess, sys
from pathlib import Path
from plasmid_evidence import fasta
import finder_adapter as finder
VERSION = '1.0.0'

def union_length(regions):
    total=0; last=0
    for lo,hi in sorted(regions):
        total+=max(0,hi-max(lo-1,last));last=max(last,hi)
    return total

def circular_segments(graph, records):
    links=set(); neighbors={}; segments={}
    for line in graph.read_text().splitlines():
        f=line.split('\t')
        if f[0]=='S': segments[f[1]]=f[2].upper()
        elif f[0]=='L' and len(f)>=6:
            if f[1]==f[3] and f[2]==f[4]:links.add(f[1])
            elif f[1]!=f[3]:
                neighbors.setdefault(f[1],set()).add(f[3]);neighbors.setdefault(f[3],set()).add(f[1])
    for n,s in records.items():
        if n in segments and segments[n]!='*' and segments[n]!=s:raise ValueError('Graph DNA mismatch: '+n)
    return {n for n in links if n in records and not neighbors.get(n)}

def cigar_crosses(cigar, start, join, flank):
    blocks=re.findall(r'(\d+)([MIDNSHP=X])',cigar)
    if not blocks or ''.join(n+op for n,op in blocks)!=cigar: raise ValueError('Malformed CIGAR')
    pos=start
    for size,op in blocks:
        size=int(size)
        if op in 'M=X':
            end=pos+size
            if pos<=join-flank and end>=join+flank:return True
            pos=end
        elif op in 'DN':pos+=size
    return False

def junction_counts(lines, reference, join, flank, minmq):
    reads=0; perfect=0; names=set()
    for line in lines:
        f=line.rstrip().split('\t')
        if len(f)<11:raise ValueError('Malformed SAM row')
        flag=int(f[1])
        if flag & 3844 or int(f[4])<minmq or f[2]!=reference:continue
        if cigar_crosses(f[5],int(f[3])-1,join,flank):
            reads+=1;names.add(f[0])
            if 'NM:i:0' in f[11:] and not re.search('[IDNSH]',f[5]):perfect+=1
    return dict(crossing_reads=reads,distinct_pair_names=len(names),perfect_unclipped_crossing_reads=perfect)

def candidate_matches(path, candidates, assemblies, minidentity, mincoverage):
    grouped={}
    for line in path.read_text().splitlines():
        if not line.strip():continue
        f=line.split('\t')
        if len(f)!=8:raise ValueError('Unexpected local BLAST columns')
        q,s=f[:2]
        if q not in candidates or s not in assemblies:raise ValueError('Local BLAST sequence ID mismatch')
        identity=float(f[2]);lo,hi=sorted((int(f[3]),int(f[4])))
        if int(f[5])!=len(candidates[q]) or not 1<=lo<=hi<=len(candidates[q]):raise ValueError('Invalid local BLAST coordinates')
        if identity>=minidentity:grouped.setdefault((q,s),[]).append((lo,hi))
    rows=[]
    for (q,s),regions in grouped.items():
        cov=100*union_length(regions)/len(candidates[q])
        if cov>=mincoverage:rows.append(dict(candidate=q,reconstructed_id=s,query_coverage_percent=cov))
    return rows

def run(cmd,log,commands,stdout=None):
    commands.append(cmd);print('Running: '+' '.join(cmd),flush=True)
    with log.open('w') as err:
        if stdout:
            with stdout.open('w') as out:subprocess.run(cmd,stdout=out,stderr=err,check=True)
        else:subprocess.run(cmd,stdout=err,stderr=err,check=True)

def spades_launcher(path, executable, memory):
    # Unicycler 0.5 appends -m 1024. Enforce one SPAdes memory argument.
    source='''#!/usr/bin/env python3
import os, sys
args=[]; i=1
while i<len(sys.argv):
    a=sys.argv[i]
    if a in ('-m','--memory'):
        i+=2; continue
    if a.startswith('--memory='):
        i+=1; continue
    args.append(a); i+=1
if any(a in ('-o','--output') or a.startswith('--output=') for a in args):
    args += ['-m', MEMORY]
os.execv(EXECUTABLE, [EXECUTABLE]+args)
'''.replace('MEMORY',repr(str(memory))).replace('EXECUTABLE',repr(str(executable)))
    path.write_text(source);path.chmod(0o755)

def mapping(bwa,samtools,reference,r1,r2,out,threads,commands):
    run([bwa,'index',str(reference)],out/'bwa_index.log',commands)
    bam=out/'junction.sorted.bam'
    first=[bwa,'mem','-t',str(threads),str(reference),str(r1),str(r2)]
    second=[samtools,'sort','-@','1','-m','256M','-T',str(out/'sort_tmp'),'-o',str(bam),'-']
    commands.extend([first,second]);print('Mapping reads across proposed junction',flush=True)
    with (out/'bwa.log').open('w') as l1,(out/'samtools_sort.log').open('w') as l2:
        proc=subprocess.Popen(first,stdout=subprocess.PIPE,stderr=l1)
        try:
            sorted_proc=subprocess.run(second,stdin=proc.stdout,stderr=l2)
            proc.stdout.close();code=proc.wait()
        finally:
            if proc.poll() is None:proc.terminate();proc.wait()
        if code or sorted_proc.returncode:raise ValueError('Junction mapping or sorting failed; inspect logs')
    run([samtools,'quickcheck','-v',str(bam)],out/'quickcheck.log',commands)
    run([samtools,'index',str(bam)],out/'index.log',commands)
    return bam

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version',action='version',version=VERSION)
    p.add_argument('--assembly',required=True,type=Path)
    p.add_argument('--reads1',required=True,type=Path);p.add_argument('--reads2',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--rrna-gff',type=Path);p.add_argument('--mob-report',type=Path);p.add_argument('--gfa',type=Path)
    p.add_argument('--run-barrnap',action='store_true');p.add_argument('--run-mob',action='store_true')
    p.add_argument('--reuse-unicycler',type=Path,help='Reuse instead of reconstructing; junction mapping still runs')
    p.add_argument('--recover-from',type=Path,help='Explicit partial recovery from an interrupted workflow; uses a NEW output folder')
    p.add_argument('--unicycler-python',default='/usr/bin/python3',help='Python with installed Unicycler 0.5.0 for graph recovery')
    p.add_argument('--threads',type=int,default=2);p.add_argument('--memory-gb',type=int,default=7)
    p.add_argument('--junction-arm',type=int,default=500);p.add_argument('--min-anchor',type=int,default=20)
    p.add_argument('--min-mapq',type=int,default=30)
    p.add_argument('--match-identity',type=float,default=95);p.add_argument('--match-coverage',type=float,default=80)
    p.add_argument('--annotation-genbank',type=Path)
    p.add_argument('--annotation-gff',type=Path)
    p.add_argument('--annotation-fasta',type=Path)
    p.add_argument('--cargo-config',type=Path,help='Optional functional cargo modules; independent of plasmid calls')
    finder.arguments(p)
    a=p.parse_args()
    finder.validate(a)
    if a.threads<1 or a.memory_gb<1 or a.min_anchor<1 or a.junction_arm<a.min_anchor:p.error('Invalid resource or junction settings')
    if not 0<=a.min_mapq<=255 or not 0<a.match_identity<=100 or not 0<a.match_coverage<=100:p.error('Invalid quality or match settings')
    for path in [a.assembly,a.reads1,a.reads2,a.rrna_gff,a.mob_report,a.gfa,a.annotation_gff,a.annotation_fasta,a.annotation_genbank]:
        if path and not path.is_file():p.error('Missing file: '+str(path))
    if a.reads1.resolve()==a.reads2.resolve():p.error('Read mates must be different files')
    if a.recover_from and a.reuse_unicycler:p.error('Choose either --recover-from or --reuse-unicycler')
    if a.recover_from:
        previous=json.loads((a.recover_from/'workflow.json').read_text())
        old=previous.get('settings',{})
        for key in ['assembly','reads1','reads2']:
            if key not in old or Path(old[key]).resolve()!=getattr(a,key).resolve():p.error('Recovery input differs: '+key)
        saved=previous.get('inputs',{}).get('assembly',{}).get('sha256')
        if not saved or saved!=hashlib.sha256(a.assembly.read_bytes()).hexdigest():p.error('Recovery assembly fingerprint differs or is missing')
        for r,record in zip([a.reads1,a.reads2],previous.get('read_files',[])):
            if r.stat().st_size!=record['size_bytes'] or r.stat().st_mtime_ns!=record['mtime_ns']:p.error('Read file changed since original run')
        if len(previous.get('read_files',[]))!=2:p.error('Missing original read metadata')
    if a.out.exists():p.error('Choose a new output folder; existing results are preserved')
    a.out.mkdir(parents=True)
    commands=[];metadata={'version':VERSION,'settings':{k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},'commands':commands,'state':'started'}
    def checkpoint(): (a.out/'workflow.json').write_text(json.dumps(metadata,indent=2)+'\n')
    checkpoint()
    try:
        inputs={}
        for key in ['assembly','rrna_gff','mob_report','gfa']:
            path=getattr(a,key)
            if path:inputs[key]={'path':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
        metadata['inputs']=inputs
        metadata['read_files']=[{'path':str(r.resolve()),'size_bytes':r.stat().st_size,'mtime_ns':r.stat().st_mtime_ns} for r in [a.reads1,a.reads2]]
        screen=[sys.executable,str(Path(__file__).with_name('plasmid_evidence.py')),'--assembly',str(a.assembly.resolve()),'--out',str((a.out/'screen').resolve()),'--threads',str(a.threads)]
        for arg in ['rrna_gff','mob_report','gfa','annotation_gff','annotation_fasta','annotation_genbank']:
            if getattr(a,arg):screen += ['--'+arg.replace('_','-'),str(getattr(a,arg).resolve())]
        for arg in ['run_barrnap','run_mob']:
            if getattr(a,arg):screen+=['--'+arg.replace('_','-')]
        for arg in ['plasmidfinder_results','plasmidfinder_db']:
            if getattr(a,arg): screen += ['--'+arg.replace('_','-'),str(getattr(a,arg).resolve())]
        for arg in ['plasmidfinder_executable','mob_executable','plasmidfinder_identity','plasmidfinder_coverage']:
            screen += ['--'+arg.replace('_','-'),str(getattr(a,arg))]
        if a.plasmidfinder_python: screen += ['--plasmidfinder-python',a.plasmidfinder_python]
        if a.run_plasmidfinder: screen.append('--run-plasmidfinder')
        run(screen,a.out/'screen.log',commands)
        metadata['annotation_evidence']=json.loads((a.out/'screen/annotation_evidence.json').read_text())
        metadata['finder_hits']=json.loads((a.out/'screen/plasmidfinder_hits.json').read_text())
        metadata['screen_provenance']=json.loads((a.out/'screen/provenance.json').read_text())
        candidatefile=a.out/'screen/candidates.fasta'
        if not candidatefile.stat().st_size:
            metadata['state']='no_candidates';checkpoint()
            (a.out/'report.md').write_text('# Reconstruction workflow\n\nNo candidates in supplied screening checks. Reconstruction was skipped; plasmid absence is not established.\n')
            print('No candidates; reconstruction skipped.');return
        executables={}
        required=['bwa','samtools','blastn','makeblastdb']+([] if a.reuse_unicycler else ['unicycler','spades'])
        for name in required:
            resolved=shutil.which(name)
            if not resolved:raise ValueError('Missing executable: '+name)
            executables[name]=str(Path(resolved).resolve())
        metadata['executables']=executables
        for name,exe in executables.items():
            versionargs=['--version'] if name in ['unicycler','spades','samtools'] else ['-version'] if name in ['blastn','makeblastdb'] else []
            info=subprocess.run([exe]+versionargs,capture_output=True)
            (a.out/(name+'_version.txt')).write_bytes(info.stdout+info.stderr)
        u=a.reuse_unicycler.resolve() if a.reuse_unicycler else (a.out/'unicycler').resolve()
        if not a.reuse_unicycler:
            launcher=(a.out/'spades_memory_launcher.py').resolve()
            spades_launcher(launcher,executables['spades'],a.memory_gb)
            if a.recover_from:
                source=(a.recover_from/'unicycler').resolve()
                run([a.unicycler_python,str(Path(__file__).with_name('recover_graphs.py')),'--source',str(source),'--out',str(u),'--reads1',str(a.reads1.resolve()),'--reads2',str(a.reads2.resolve())],a.out/'graph_recovery.log',commands)
                metadata['recovery']=json.loads((u/'recovery_notes.json').read_text());checkpoint()
            run([executables['unicycler'],'-1',str(a.reads1.resolve()),'-2',str(a.reads2.resolve()),'-o',str(u),'-t',str(a.threads),'--mode','normal','--spades_path',str(launcher),'--keep','3'],a.out/'unicycler_terminal.log',commands)
        assemblies=fasta(u/'assembly.fasta');candidates=fasta(candidatefile)
        circles=circular_segments(u/'assembly.gfa',assemblies)
        metadata['state']='reconstructed';metadata['unicycler_directory']=str(u);checkpoint()
        db=a.out/'reconstructed_db'
        run([executables['makeblastdb'],'-in',str(u/'assembly.fasta'),'-dbtype','nucl','-out',str(db)],a.out/'makeblastdb.log',commands)
        matchesfile=a.out/'candidate_matches.tsv'
        run([executables['blastn'],'-query',str(candidatefile),'-db',str(db),'-task','blastn','-dust','no','-max_target_seqs',str(len(assemblies)),'-outfmt','6 qseqid sseqid pident qstart qend qlen length bitscore','-out',str(matchesfile)],a.out/'local_blast.log',commands)
        matches=candidate_matches(matchesfile,candidates,assemblies,a.match_identity,a.match_coverage)
        metadata['candidate_matches']=matches;metadata['isolated_circular_ids']=sorted(circles)
        targets=sorted({m['reconstructed_id'] for m in matches if m['reconstructed_id'] in circles})
        results=[]
        for n in targets:
            safe=hashlib.sha256(n.encode()).hexdigest()[:12];folder=a.out/('junction_'+safe);folder.mkdir()
            sequence=assemblies[n];arm=min(a.junction_arm,len(sequence)//2)
            if arm<a.min_anchor:
                results.append({'reconstructed_id':n,'status':'too short for specified junction anchors'});continue
            refname='junction_'+safe;junction=folder/'junction.fasta'
            junction.write_text('>'+refname+'\n'+sequence[-arm:]+sequence[:arm]+'\n')
            if fasta(junction)[refname]!=sequence[-arm:]+sequence[:arm]:raise ValueError('Junction verification failed')
            reference=junction
            bam=mapping(executables['bwa'],executables['samtools'],reference,a.reads1.resolve(),a.reads2.resolve(),folder,a.threads,commands)
            cmd=[executables['samtools'],'view','-q',str(a.min_mapq),'-F','3844',str(bam),f'{refname}:{arm-a.min_anchor+1}-{arm+a.min_anchor}']
            commands.append(cmd)
            with (folder/'count.log').open('w') as log:
                proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=log,text=True)
                try:counts=junction_counts(proc.stdout,refname,arm,a.min_anchor,a.min_mapq)
                finally:
                    proc.stdout.close()
                    if proc.poll() is None:code=proc.wait()
                    else:code=proc.returncode
                if code:raise ValueError('samtools junction count failed')
            counts['mapping_reference']='junction only'
            counts.update(reconstructed_id=n,length_bp=len(sequence),junction_arm_bp=arm,min_anchor_bp=a.min_anchor,min_mapq=a.min_mapq,excluded_flags=3844,
                          status='circular candidate with junction-read support' if counts['crossing_reads'] else 'circular graph candidate; no qualifying junction reads',
                          limitation='Junction-only mapping does not test genome-wide uniqueness, PCR independence or plasmid identity.')
            (folder/'support.json').write_text(json.dumps(counts,indent=2)+'\n');results.append(counts)
            if counts['crossing_reads']:
                (folder/'supported_circular_candidate.fasta').write_text('>'+n+'\n'+sequence+'\n')
        if a.cargo_config:
            selected={n:assemblies[n] for n in targets}
            circular_sources={m['candidate'] for m in matches if m['reconstructed_id'] in targets}
            for n,sequence in candidates.items():
                if n not in circular_sources:
                    if n in selected and selected[n]!=sequence:raise ValueError('Candidate cargo ID collision')
                    selected[n]=sequence
            cargo_fasta=a.out/'cargo_candidates.fasta'
            cargo_fasta.write_text(''.join('>'+n+'\n'+seq+'\n' for n,seq in selected.items()))
            metadata['cargo_sequence_status']={n:'isolated circular graph candidate' if n in targets else 'unresolved original candidate' for n in selected}
            from cargo_analysis import analyze
            metadata['cargo_analysis']=analyze(cargo_fasta,a.cargo_config,a.out/'cargo')
        metadata['results']=results;metadata['state']='completed';checkpoint()
        report=['# Active reconstruction v'+VERSION,'Original review candidates: '+str(len(candidates)),
                'Reconstruction: '+('partial saved-graph recovery (not exact resume)' if a.recover_from else 'existing run reused' if a.reuse_unicycler else 'Unicycler run from paired reads'),
                'Circular matches: '+str(len(targets)),
                'Matching cutoffs are exploratory tracking rules, not validated plasmid classification. Multi-segment cycles are not assessed. Unmatched candidates remain unresolved.',
                'Junction reads are mapped to a synthetic end-to-start reference only. Counts are supporting evidence, not proof of plasmid identity.']
        if a.recover_from:report.append('Recovery details: '+json.dumps(metadata['recovery'],indent=2))
        report.append('## Screening evidence (original contig IDs)')
        with (a.out/'screen/evidence.tsv').open() as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                if row['status'] != 'no positive evidence in supplied checks':
                    report.append(row['contig']+': '+row['status']+'. '+row['explanation'])
        report += [json.dumps(r,indent=2) for r in results]
        if not targets:report.append('No isolated circular component met candidate matching cutoffs. This does not establish plasmid absence.')
        (a.out/'report.md').write_text('\n\n'.join(report)+'\n');print('Completed. Read',a.out/'report.md')
    except BaseException as error:
        metadata['state']='interrupted' if isinstance(error,KeyboardInterrupt) else 'failed';metadata['error']=str(error);checkpoint();raise

if __name__=='__main__':
    try:main()
    except KeyboardInterrupt:print('Interrupted; logs and partial outputs retained.',file=sys.stderr);sys.exit(130)
    except (ValueError,OSError,subprocess.CalledProcessError) as error:print('ERROR:',error,file=sys.stderr);sys.exit(1)
