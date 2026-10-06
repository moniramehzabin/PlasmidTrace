# PlasmidTrace

Experimental prototype **v1.0.0** for screening, reconstructing and investigating plasmid candidates from bacterial short-read assemblies.

PlasmidTrace follows candidate sequences from the initial assembly through reconstruction, checks read support across proposed circular junctions, and optionally investigates their functional genes.

Biological benchmarking is pending. Results are evidence for review, not confirmed plasmid classifications.

## Workflow

1. **Screen candidates:** examine assembly-header coverage, predicted rRNA overlaps, optional PlasmidFinder and MOB-suite results, annotation clues and matching-ID graph self-links.
2. **Reconstruct:** run Unicycler from paired reads, reuse a completed assembly, or recover completed saved graphs from an interrupted run.
3. **Track candidates:** use local BLAST to identify where the original candidates occur in the reconstructed assembly.
4. **Check circularity:** identify isolated single-segment circular graph candidates and count reads crossing a synthetic end-to-start junction.
5. **Analyze functional cargo:** optionally run or import AMRFinderPlus, VirulenceFinder and Bioremediation Gene Miner results for selected candidate sequences.

## Requirements

- Python 3.10 or newer.
- External programs and databases installed separately.
- Screening, depending on selected options: barrnap, MOB-suite and PlasmidFinder.
- Reconstruction: Unicycler 0.5.0, SPAdes, BLAST, BWA and SAMtools.
- Functional analysis: AMRFinderPlus, VirulenceFinder and/or Bioremediation Gene Miner with their required databases and environments.

The PlasmidTrace Python scripts use the standard library. External modules may require their own Python dependencies.

Genomic datasets and external databases are not included.

## Getting started

Download this repository, open a terminal in its directory, and choose the appropriate command below.

### Screen an assembly

```bash
python3 plasmid_evidence.py \
  --assembly assembly.fasta \
  --run-barrnap \
  --threads 2 \
  --out screen_results
```

This screens the assembly and runs barrnap to assess rRNA overlaps. Additional evidence sources can be supplied through optional arguments.

### Reconstruct and check candidates

```bash
python3 reconstruct.py \
  --assembly assembly.fasta \
  --reads1 R1.fastq.gz \
  --reads2 R2.fastq.gz \
  --run-barrnap \
  --threads 2 \
  --memory-gb 7 \
  --out reconstruction_results
```

This screens candidates, runs Unicycler from the paired reads, tracks candidates into the reconstruction and checks eligible circular junctions.

Reconstruction runs on the supplied paired reads and can take several hours. Finding a candidate does not guarantee that reconstruction will resolve it into a circle.

### Analyze functional cargo

```bash
python3 cargo_analysis.py \
  --assembly candidate.fasta \
  --config cargo_config.json \
  --out cargo_results
```

Use a FASTA containing the candidate sequences to investigate. Configuration can specify existing result files or commands that run external modules.

See [CARGO_MODULES.md](CARGO_MODULES.md), [cargo_import.example.json](cargo_import.example.json) and [cargo_run.example.json](cargo_run.example.json).

Use a new or empty output directory for each run. Commands configured for external modules are executed without a shell; configure only trusted programs.

For all available options:

```bash
python3 plasmid_evidence.py --help
python3 reconstruct.py --help
python3 cargo_analysis.py --help
```

## Candidate screening

PlasmidTrace retains the evidence supporting each candidate rather than treating every positive clue as confirmation.

- **Coverage:** elevated coverage relative to an exploratory assembly baseline can flag a candidate. This is original SPAdes header coverage, not mapped base depth or plasmid copy number.
- **rRNA:** substantial overlap with predicted rRNA produces an rRNA-dominated warning. This takes priority over a positive candidate call.
- **Plasmid markers:** optional PlasmidFinder and MOB-suite evidence can contribute to screening.
- **Annotation:** supported replication or partition-related annotation can contribute evidence; annotation alone does not prove plasmid identity.
- **Graph structure:** an isolated same-orientation self-link is distinguished from a self-link connected to other segments.

Missing checks remain unassessed. Absence of candidates does not establish plasmid absence.

## Reconstruction and recovery

The reconstruction workflow can:

- Run Unicycler from paired reads.
- Reuse an existing completed Unicycler assembly.
- Recover and score completed saved graphs from an interrupted run.

Partial graph recovery is **not an exact SPAdes resume**. It uses completed saved graphs and does not finish missing k-mer stages. Recovery details are recorded in the output.

Candidate-to-reconstruction matches use exploratory identity and query-coverage cutoffs. These track sequences; they are not validated plasmid classification thresholds.

Multi-segment circular paths are not assessed by this prototype.

## Junction-read support

For an eligible circular candidate, PlasmidTrace constructs a reference joining the sequence end to its beginning, maps reads to that reference, and counts qualifying alignments crossing the join.

The output reports:

- Crossing reads.
- Distinct read-pair names.
- Perfect, unclipped crossing reads.

Junction support strengthens evidence for the proposed closure but does not confirm plasmid identity.

Mapping quality is relative to a junction-only reference. Genome-wide uniqueness is not established, and distinct read-pair names do not establish independent original DNA molecules. Read counts are reported as support without a validated confidence threshold.

## Candidate-only functional analysis

Adding `--cargo-config` to `reconstruct.py` runs functional analysis after reconstruction and junction checks.

The generated `cargo_candidates.fasta` contains:

- Isolated circular candidates matched to screened candidates.
- Original screened candidates that remain unresolved.

It does not automatically include the whole assembly. If no candidates are screened, cargo analysis is skipped.

Candidate statuses remain in `workflow.json`. Functional results do not change plasmid status or prove plasmid carriage.

### AMRFinderPlus

Reports resistance-associated and other supported elements. Element types are preserved separately, so stress and virulence entries are not automatically counted as AMR genes.

### VirulenceFinder

Reports matches to the selected virulence databases. Database coverage depends on taxon; no hit does not establish absence of virulence.

### Bioremediation Gene Miner

Investigates candidate bioremediation-associated functions while preserving confidence labels and review or reject decisions.

Gene Miner requires matching GenBank annotation. Its annotated nucleotide sequences must match the selected candidate DNA exactly or by reverse complement.

PlasmidTrace does not automatically annotate newly reconstructed circles. Whole-genome annotation cannot be used directly against a candidate-only FASTA unless the required matching candidate annotation is supplied.

An empty Gene Miner TSV requires an explicit workbook path and verification of an empty `All_candidates` sheet; otherwise it is treated as an error. Original TSV bytes and the workbook hash are retained.

Predicted genes do not establish resistance, pathogenicity or bioremediation activity.

## Outputs and provenance

Depending on the workflow and selected modules, outputs include:

- Candidate FASTA files and evidence tables.
- Human-readable reports.
- Candidate-to-reconstruction matches.
- Junction-read support results.
- Combined functional cargo tables.
- Raw external-tool results and execution logs.
- Workflow settings and provenance records.

Provenance includes recorded commands, input information and applicable file hashes. Retain these outputs when reporting an analysis.

## Validation

Run the software tests with:

```bash
python3 -m unittest discover -q
```

The current release passes **40 software tests**, covering parsing, coordinate checks, candidate tracking, graph context, recovery, empty-result handling and simulated workflow execution.

The GitHub Actions test workflow has also passed.

Real candidate-only tests have exercised zero-result imports and automatic execution of all three cargo modules.

The **AMRFinderPlus positive control reported 15 hits**. All 15 were imported with contig IDs, gene names, element types, methods and coordinates preserved.

Positive-hit controls for the other cargo modules and end-to-end candidate-only cargo execution remain to be verified.

Software tests and adapter checks do not measure biological sensitivity or specificity. Biological benchmarking and comparisons with other plasmid tools are pending. No claim of uniqueness or superior performance is made.

## Interpretation limits

- Elevated coverage alone does not confirm a plasmid.
- Graph circularity and junction-read support do not confirm plasmid identity.
- Marker and reference databases have limited coverage.
- Chromosome reference matches do not automatically exclude plasmids.
- Functional cargo results do not establish plasmid carriage.
- No detected hits do not prove absence of plasmids, resistance, virulence or bioremediation functions.

## External software

PlasmidTrace uses or supports evidence from:

- [Unicycler](https://github.com/rrwick/Unicycler)
- [MOB-suite](https://github.com/phac-nml/mob-suite)
- [PlasmidFinder](https://github.com/genomicepidemiology/plasmidfinder)
- [AMRFinderPlus](https://github.com/ncbi/amr)
- [VirulenceFinder](https://github.com/genomicepidemiology/virulencefinder)
- [Bioremediation Gene Miner](https://github.com/moniramehzabin/Bioremediation-Gene-Miner)

External software and databases retain their own licenses and citation requirements. Cite the external programs and databases used in your analysis.

## Authors and citation

**Monira Mehzabin and Khandoker Md Rezwan**

Citation metadata is provided in [CITATION.cff](CITATION.cff). GitHub also provides a **Cite this repository** option.

## License

PlasmidTrace is released under the [MIT license](LICENSE).

External tools and databases are not covered by the PlasmidTrace license.
