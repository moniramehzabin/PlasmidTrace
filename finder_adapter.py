"""Strict adapter for PlasmidFinder results_tab.tsv (assembly BLAST mode)."""
import csv
import hashlib
import re
import shutil
from pathlib import Path

def arguments(parser):
    parser.add_argument('--plasmidfinder-python',help='Python executable containing the plasmidfinder module (PF3)')
    parser.add_argument('--finders', action='store_true', help='Run MOB-Recon and PlasmidFinder; requires configured installations and database')
    parser.add_argument('--run-plasmidfinder', action='store_true')
    parser.add_argument('--plasmidfinder-results', type=Path, help='Import results_tab.tsv from a completed assembly search')
    parser.add_argument('--plasmidfinder-db', type=Path)
    parser.add_argument('--plasmidfinder-executable', default='plasmidfinder.py')
    parser.add_argument('--mob-executable', default='mob_recon')
    parser.add_argument('--plasmidfinder-identity', type=float, default=.95)
    parser.add_argument('--plasmidfinder-coverage', type=float, default=.60)

def validate(a):
    if a.finders: a.run_mob = a.run_plasmidfinder = True
    if a.run_plasmidfinder and a.plasmidfinder_results:
        raise ValueError('Use --run-plasmidfinder OR --plasmidfinder-results')
    if not 0 < a.plasmidfinder_identity <= 1 or not 0 < a.plasmidfinder_coverage <= 1:
        raise ValueError('PlasmidFinder thresholds must be fractions in (0,1]')
    if a.run_plasmidfinder:
        if not a.plasmidfinder_db or not a.plasmidfinder_db.is_dir():
            raise ValueError('Provide installed database directory with --plasmidfinder-db')
        if not shutil.which(a.plasmidfinder_python or a.plasmidfinder_executable):
            raise ValueError('Missing PlasmidFinder executable: '+a.plasmidfinder_executable)
    if a.run_mob and not shutil.which(a.mob_executable):
        raise ValueError('Missing MOB-Recon executable: '+a.mob_executable)

def database_fingerprint(path):
    if not path: return {}
    return {str(p.relative_to(path)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(path.rglob('*')) if p.is_file() and p.suffix in ('.fsa','.fasta','.fa')}

def import_hits(path, records, identity, coverage):
    hits = {n: [] for n in records}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'Plasmid','Identity','Query / Template length','Contig'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError('Unsupported PlasmidFinder table: expected '+', '.join(sorted(required)))
        for row in reader:
            if not any(row.values()): continue
            name = row['Contig'].strip()
            if name not in records: raise ValueError('PlasmidFinder contig ID mismatch: '+name)
            pct = float(row['Identity'].rstrip('%'))
            m = re.fullmatch(r'\s*(\d+)\s*/\s*(\d+)\s*', row['Query / Template length'])
            if not m: raise ValueError('Invalid PlasmidFinder marker length ratio')
            aligned, template = map(int,m.groups())
            if template <= 0 or aligned <= 0 or not 0 <= pct <= 100:
                raise ValueError('Invalid PlasmidFinder hit values')
            cov = min(100.0,100*aligned/template)
            if pct >= 100*identity and cov >= 100*coverage:
                hits[name].append(dict(marker=row['Plasmid'],identity_percent=pct,
                    marker_coverage_percent=cov,accession=row.get('Accession number',row.get('Accession','')),
                    position=row.get('Position in contig','')))
    return hits

def command(a, assembly, out):
    prefix = [a.plasmidfinder_python,'-m','plasmidfinder'] if a.plasmidfinder_python else [a.plasmidfinder_executable]
    return prefix+['-x','-i',str(assembly.resolve()),'-o',str(out.resolve()),
            '-p',str(a.plasmidfinder_db.resolve()),'-t',str(a.plasmidfinder_identity),
            '-l',str(a.plasmidfinder_coverage)]
