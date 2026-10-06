#!/usr/bin/env python3
"""Explicit partial recovery of completed Unicycler 0.5.0 SPAdes graphs."""
import argparse, hashlib, json, re, shutil, statistics, subprocess
from pathlib import Path

def inventory(source):
    rows=[]
    for path in sorted(source.glob('001_spades_graph_k*.gfa')):
        match=re.fullmatch(r'001_spades_graph_k(\d+)\.gfa',path.name)
        if match and path.stat().st_size:
            rows.append((int(match[1]),path))
    if not rows: raise ValueError('No completed saved SPAdes graphs found')
    return rows

def insert_sizes(log):
    values=re.findall(r'Insert size = ([\d.]+), deviation = ([\d.]+)',log)
    if not values: raise ValueError('No saved insert-size estimates; refusing guessed metadata')
    return tuple(statistics.median(float(x[i]) for x in values) for i in (0,1))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--reads1',required=True);p.add_argument('--reads2',required=True)
    a=p.parse_args()
    import unicycler.spades_func as sf
    installed=subprocess.run(['unicycler','--version'],capture_output=True,text=True,check=True).stdout.strip()
    if installed not in ('Unicycler v0.5.0','Unicycler version: v0.5.0'): raise ValueError('Graph recovery supports installed Unicycler 0.5.0 only')
    rows=inventory(a.source)
    logs=[a.source/'spades_assembly/spades.log',a.source.parent/'shutdown_checkpoint/spades.log']
    # Prefer preserved pre-interruption metadata if available.
    if logs[1].is_file(): logs.reverse()
    estimates=None
    for log in logs:
        if log.is_file():
            try: estimates=insert_sizes(log.read_text(errors='replace'));break
            except ValueError: pass
    if estimates is None: raise ValueError('No usable insert-size log found')
    if a.out.exists(): raise ValueError('Recovery output already exists; choose a new folder')
    a.out.mkdir(parents=True)
    copies=[]; hashes={}
    for k,path in rows:
        dest=a.out/path.name;shutil.copy2(path,dest);copies.append(str(dest))
        hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    kmers=[k for k,_ in rows]
    sf.get_kmer_range=lambda *args,**kwargs: kmers
    sf.run_spades_all_kmers=lambda *args,**kwargs: (copies,*estimates)
    graph=sf.get_best_spades_graph(a.reads1,a.reads2,None,str(a.out),0.25,1,
        '/usr/bin/spades',2,3,len(kmers),0.2,0.95,kmers,0,False,str(a.out/'001_spades_graph'),None)
    graph.save_to_gfa(str(a.out/'002_depth_filter.gfa'),save_copy_depth_info=True,newline=True,include_insert_size=True)
    notes={'mode':'partial graph recovery; not exact SPAdes resume','completed_kmers_used':kmers,
        'source':str(a.source.resolve()),'source_graph_sha256':hashes,'insert_size_estimates':estimates,
        'limitation':'Only completed saved graphs were scored; unfinished k-mer stages were not run.'}
    (a.out/'recovery_notes.json').write_text(json.dumps(notes,indent=2)+'\n')
    print('Prepared partial graph recovery:',a.out)

if __name__=='__main__': main()
