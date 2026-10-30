"""Describe checks actually assessed; a successful process is not biological proof."""
import csv,json
from pathlib import Path

def summarize(folder,settings):
    folder=Path(folder);path=folder/'output/evidence.tsv'
    if not path.exists():path=folder/'output/screen/evidence.tsv'
    if not path.exists():
        return {'limited':False,'label':'Workflow completed','warnings':['Review workflow reports; completion does not establish plasmid identity.'],'checks':{}}
    with path.open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    def assessed(field):return bool(rows) and any(r.get(field,'not assessed')!='not assessed' for r in rows)
    checks={label:assessed(field) for label,field in [('coverage','coverage_ratio'),('rRNA','rrna_fraction'),('MOB-Recon','mob_classification'),('PlasmidFinder','plasmidfinder_markers'),('graph','graph_context'),('annotation','annotation_gene_evidence')]}
    warnings=[]
    if not checks['MOB-Recon'] and not checks['PlasmidFinder']:warnings.append('No plasmid-marker checks were assessed. Zero candidates does not mean plasmids are absent.')
    if not checks['coverage']:warnings.append('Coverage could not be assessed from these FASTA headers. Original SPAdes coverage values were not supplied.')
    if not checks['rRNA']:warnings.append('rRNA evidence was not assessed.')
    warnings.append('Computational evidence supports review, not confirmed plasmid identity or biological absence.')
    limited=not checks['coverage'] or not checks['rRNA'] or not (checks['MOB-Recon'] and checks['PlasmidFinder'])
    return {'limited':limited,'label':'Run completed · limitations noted' if limited else 'Screening checks completed','checks':checks,'contigs':len(rows),'review_candidates':sum(r.get('status')=='candidate for review' for r in rows),'warnings':warnings}
