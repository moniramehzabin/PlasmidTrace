# PlasmidTrace

Experimental prototype v1.0.0 for screening, reconstructing and investigating
plasmid candidates from bacterial short-read assemblies.

## Workflow

1. Screen coverage, rRNA overlaps, optional PlasmidFinder/MOB-suite results,
   annotation clues and matching-ID graph self-links.
2. Run Unicycler from paired reads, reuse a completed run, or explicitly recover
   completed saved graphs from an interrupted run.
3. Track original candidates into reconstruction by local BLAST.
4. Identify isolated single-segment circular graph candidates and count reads
   crossing a synthetic end-to-start junction.
5. Optionally run/import AMRFinderPlus, VirulenceFinder and Bioremediation Gene
   Miner results on candidate sequences.

Coverage alone is exploratory. Graph circularity and junction support do not
confirm plasmid identity. Marker/reference databases have limited coverage.
No claim of uniqueness or superior performance is made; benchmarking is pending.
Competitive mapping and alternative-junction signature auditing are excluded.

## Installation

Python 3.10+; the adapters use Python's standard library. External programs and
databases must be installed separately. Screening may need barrnap, MOB-suite,
PlasmidFinder. Reconstruction needs Unicycler 0.5.0, SPAdes, BLAST, BWA and SAMtools.
Gene Miner requires its own environment, annotation and curated database.
No genomic datasets or external databases are included.

## Examples

    python3 plasmid_evidence.py --assembly assembly.fasta --run-barrnap --threads 2 --out screen_results

    python3 reconstruct.py --assembly assembly.fasta --reads1 R1.fastq.gz --reads2 R2.fastq.gz --run-barrnap --threads 2 --memory-gb 7 --out reconstruction_results

    python3 cargo_analysis.py --assembly candidate.fasta --config cargo_config.json --out cargo_results

Use new output folders. Read CARGO_MODULES.md for gene-analysis configuration.
Explicit argv commands are run without a shell. Only configure trusted programs.

## Candidate-only cargo analysis

Adding --cargo-config to reconstruct.py runs cargo analysis after reconstruction
and junction checks. cargo_candidates.fasta contains isolated circular candidates
matched to screened candidates, plus original screened candidates that remain
unresolved. It never contains the whole assembly. workflow.json retains the
candidate statuses; cargo tables by themselves do not establish plasmid carriage.
If no candidates are screened, cargo analysis is skipped.

GenBank annotation supplied to Gene Miner must match the selected candidate DNA
exactly or by reverse complement. The tool does not automatically annotate newly
reconstructed circles: prepare matching candidate annotation or run AMRFinder and
VirulenceFinder first. Original whole-genome Gene Miner annotation is not valid
against a candidate-only FASTA. Gene Miner review/reject decisions are preserved.
An empty TSV requires an explicit workbook path and an empty All_candidates sheet;
otherwise it is an error. Original TSV bytes and workbook hash are retained.

## Validation and limits

    python3 -m unittest discover -q

40 tests pass in this release, including mocked workflow execution, parsing,
coordinate checks, candidate tracking and empty-workbook verification. Mocked
executables do not measure biological accuracy. Real candidate-only tests have
exercised zero-result imports and automatic execution of all three cargo modules.
Positive-hit external-tool controls and end-to-end candidate-only cargo execution
remain to be verified. Multi-segment circular paths are not reconstructed by this
prototype. Partial graph recovery does not complete missing k-mer stages.

VirulenceFinder database taxon coverage is limited. Absence of hits never proves
absence of virulence, resistance, bioremediation activity or plasmids.

## External software

- https://github.com/rrwick/Unicycler
- https://github.com/phac-nml/mob-suite
- https://github.com/genomicepidemiology/plasmidfinder
- https://github.com/ncbi/amr
- https://github.com/genomicepidemiology/virulencefinder
- https://github.com/moniramehzabin/Bioremediation-Gene-Miner

External software and databases retain their own licenses and citation requirements.

## Authors and license

Monira Mehzabin and Khandoker Md Rezwan. This package is released under the
MIT license; see LICENSE. Citation metadata is in CITATION.cff. External tools
and databases must be cited separately and are not covered by this license.
