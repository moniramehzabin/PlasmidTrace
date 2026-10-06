"""Transparent keyword evidence from PGAP/Prokka GFF3; no sequence classifier."""
import re
from urllib.parse import unquote

RULES = [
    ('plasmid_replication', r'\bplasmid\b.*\breplication\b|\breplication\b.*\bplasmid\b|\breplication initiation protein\b'),
    ('plasmid_partition', r'\bplasmid\b.*\bpartition\w*\b|\bpartition\w*\b.*\bplasmid\b'),
    ('partition_context', r'\bpartition\w*\b|\bpar[abm]\b'),
    ('mobility_context', r'\brelaxase\b|\bconjugative\b|\bconjugation\b|\bmobilization\b|\bmobilisation\b'),
]

def read_gff(path, records, aliases=None):
    aliases = aliases or {}
    hits = {n: [] for n in records}
    declarations = {n: [] for n in records}
    for line in path.read_text().splitlines():
        if line.startswith('##FASTA'): break
        if not line.strip() or line.startswith('#'): continue
        f = line.split('\t')
        if len(f) != 9: raise ValueError('Invalid annotation GFF3 row')
        name = aliases.get(unquote(f[0]),unquote(f[0]))
        if name not in records: raise ValueError('Annotation ID mismatch: '+name+'; supply --annotation-fasta for exact DNA mapping')
        lo,hi = int(f[3]),int(f[4])
        if not 1 <= lo <= hi <= len(records[name]): raise ValueError('Annotation coordinates outside assembly: '+name)
        attrs = {}
        for pair in f[8].split(';'):
            if '=' in pair:
                key,value = pair.split('=',1);attrs[key] = unquote(value)
        if f[2] in ('region','source') and 'plasmid' in attrs:
            declarations[name].append(attrs['plasmid'])
        if f[2] != 'CDS': continue
        product = attrs.get('product','')
        text = product+' '+attrs.get('gene','')
        categories = [category for category,pattern in RULES if re.search(pattern,text,re.I)]
        if re.search(r'\bchromosomal replication\b|\bdnaa\b',text,re.I): categories = []
        if categories:
            hits[name].append(dict(locus_tag=attrs.get('locus_tag',attrs.get('ID','')),
                product=product,gene=attrs.get('gene',''),start=lo,end=hi,
                categories=categories,pseudo=attrs.get('pseudo','false')))
    return hits,declarations

def nominates(hits):
    return any(h['pseudo'].lower() in ('false','0','no') and
               set(h['categories']) & {'plasmid_replication','plasmid_partition'} for h in hits)

def read_genbank(path, records):
    """Read nucleotide GenBank records; refuse absent DNA or ambiguous mapping."""
    hits = {n: [] for n in records}; declarations = {n: [] for n in records}
    aliases = {}
    complement = str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')
    count = 0
    for block in path.read_text().split('\n//'):
        if not block.strip(): continue
        locus = re.search(r'^LOCUS\s+(\S+)',block,re.M)
        origin = re.search(r'^ORIGIN\s*\n([\s\S]*)',block,re.M)
        if not locus or not origin: raise ValueError('GenBank requires complete records with ORIGIN DNA')
        dna = re.sub(r'[\s\d]','',origin[1]).upper()
        if not dna or set(dna)-set('ACGTRYSWKMBDHVN'): raise ValueError('Invalid GenBank ORIGIN')
        names = [n for n,s in records.items() if dna == s or dna.translate(complement)[::-1] == s]
        if len(names) != 1: raise ValueError('GenBank DNA must match exactly one assembly record: '+locus[1])
        name = names[0]; aliases[locus[1]] = name; count += 1
        if locus[1] in records and locus[1] != name: raise ValueError('GenBank ID/DNA conflict')
        section = block.split('FEATURES',1)
        if len(section) != 2: continue
        features = re.split(r'(?m)^     (?=\S)',section[1].split('\nORIGIN',1)[0])[1:]
        for feature in features:
            first = feature.splitlines()[0].split(None,1)
            if len(first) != 2: raise ValueError('Invalid GenBank feature')
            kind,location = first
            coords = [int(x) for x in re.findall(r'\d+',location)]
            if not coords or min(coords)<1 or max(coords)>len(dna): raise ValueError('Unsupported GenBank feature location')
            attrs = {}
            for m in re.finditer(r'/(?P<key>\w+)(?:=(?:"(?P<quoted>(?:[^"]|"")*)"|(?P<plain>[^\n]+)))?',feature):
                value = m['quoted'] if m['quoted'] is not None else m['plain']
                attrs[m['key']] = 'true' if value is None else ' '.join(value.split()).replace('""','"')
            if kind == 'source' and 'plasmid' in attrs: declarations[name].append(attrs['plasmid'])
            if kind != 'CDS': continue
            product = attrs.get('product',''); text = product+' '+attrs.get('gene','')
            categories = [cat for cat,pattern in RULES if re.search(pattern,text,re.I)]
            if categories:
                hits[name].append(dict(locus_tag=attrs.get('locus_tag',''),product=product,
                    gene=attrs.get('gene',''),start=min(coords),end=max(coords),categories=categories,
                    pseudo=attrs.get('pseudo',attrs.get('pseudogene','false')),
                    annotation_record=locus[1],location=location))
    if not count: raise ValueError('No GenBank records')
    return hits,declarations,aliases
