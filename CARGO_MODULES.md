# Optional cargo modules — v1.0.0

These modules describe genes, independently of plasmid screening and circularity.
They do not nominate plasmids or change confidence in a plasmid call.
Run on the assembly whose contig names you want in the results. Reconstructed
contigs have different IDs: do not import original-assembly results against them.

## First test: import existing reports

Copy cargo_import.example.json to your own configuration, replace the paths,
and remove modules you are not using. Then run:

    python3 cargo_analysis.py --assembly YOUR_ASSEMBLY.fasta --config cargo_import.json --out NEW_CARGO_FOLDER

Results: cargo.tsv, cargo_hits.json, report.md and provenance.json. Each module
also retains the original TSV bytes. JSON retains every original column.
A new or empty output folder is required; originals are never overwritten by
import. Unknown IDs, missing headers and invalid coordinates cause an error,
rather than a false negative result. A valid header-only report means zero
reported matches in that search. Gene Miner requires the annotated GenBank used
for its report; its nucleotide sequences must match the assembly exactly or by
reverse complement. Ambiguous sequence matches and duplicate loci are refused.
Gene Miner loci are mapped to contigs, but feature coordinates are not exported.

## Run the programs from the tool

cargo_run.example.json gives explicit commands for all three programs. Replace
all installation, annotation and database paths first. Delete unused modules.
Commands are argument arrays, without a shell. {assembly} is the input FASTA;
{out} is a newly created per-module output folder. Programs run sequentially.
Each command's stdout/stderr is retained in command.log. Any failed command stops
the analysis. Programs and databases must already be installed; this adapter
neither installs dependencies nor downloads data or submits proteins online.

AMRFinderPlus uses nucleotide-only mode with --plus. This retains AMR, stress
and virulence element types separately. Providing proteins and an annotation
mapping enables additional AMRFinder analyses, but this example does not do so.
Use a locally installed database. Include database version files in
provenance_files to hash them; database identity is not established by a path alone.

VirulenceFinder example follows the tested 3.2.1 Python-module interface, using
-ifa for FASTA and -b for BLAST. The adapter imports CGE results_tab.tsv with
named columns and contig positions like 123..456.
Other formats fail explicitly. Database species coverage is limited: no hit in
an environmental isolate does not mean it lacks virulence-related genes.
The example uses 90% identity and 60% reference coverage; these are search settings,
not experimentally validated pathogenicity thresholds for your isolates.

Gene Miner is the existing Bioremediation-Gene-Miner program, not a rewritten
classifier. The command follows its v0.3.8 CLI and imports its full
Bioremediation_Gene_Miner_Report_all_candidates.tsv. It retains evidence class,
confidence, decisions, interpretation and any domain evidence in the raw JSON.
Its dependencies (including pandas, Biopython, DIAMOND and openpyxl), curated
DIAMOND DB, reference metadata and rules remain external. No automatic InterPro
web submission is enabled. Older output formats may require an explicit adapter.

## Attach to reconstruction

Add this option to your existing reconstruct.py command:

    --cargo-config /absolute/path/cargo_config.json

Cargo runs after reconstruction, on exported candidate sequences only. It includes
tracked isolated circles and unresolved original candidates, never the full assembly.
If screening finds no candidates, cargo is skipped. See README.md for annotation
matching limitations. To verify an empty Gene Miner TSV, supply a workbook setting
with the matching report path (supports {out}); populated workbooks are refused.

## Interpretation

Gene presence is a prediction. It does not establish expressed resistance,
pathogenicity, pollutant removal or plasmid carriage. Review/reject records remain
review/reject records. Missing modules are unassessed, not negative. No merged
biological score is calculated. These are optional integrations, not a claim of
novelty over existing plasmid pipelines.

## Validation

40 automated tests pass across this package, including cargo schema validation,
unknown-ID rejection, reverse-strand coordinates, preserved rejected Gene Miner
rows, and header-only versus invalid empty reports. These tests do not establish
compatibility with every released external-tool version or performance on real
isolates. Check real outputs before relying on the combined report.

Official interfaces consulted:
- https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus
- https://github.com/genomicepidemiology/virulencefinder
- https://github.com/moniramehzabin/Bioremediation-Gene-Miner
