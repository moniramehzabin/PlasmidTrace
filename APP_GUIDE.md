# PlasmidTrace browser app — preview 0.4.0

A separate local browser interface for the existing PlasmidTrace 1.0.0 scripts. Select files and options, then run; no analysis commands or JSON editing are required. The older grey desktop window has been replaced.

## Open on Windows

1. Extract the complete ZIP onto **F:**. Keep `app_assets` and all Python files together. Do not run inside the ZIP.
2. Double-click **Open_PlasmidTrace.bat**. A launcher window and a private browser page open. Keep the launcher open while analysis runs.
3. Click **Choose file** and select your genome assembly FASTA.
4. Click **Run analysis**, or **Set up and run** on the first use. First-time setup downloads the core tools and databases, then automatically starts the selected analysis. Keep the launcher open.
5. Open the report or download the results when analysis finishes.

**Advanced** contains other workflows, evidence imports, resources, output location and environment controls. The normal screening route needs none of these. An unsuccessful setup shows its log and does not start an incomplete analysis. Cancel setup to stop the pending analysis.

Windows Python 3.10+ is required. The native file picker also uses Tkinter, normally included with Windows Python. The main browser interface does not require Tkinter; a full file path can be entered if the picker is unavailable. No web framework or additional Python package is needed for the interface itself.

## Where files go

The output defaults to `results` next to the app. Extracting onto F: therefore keeps new settings, reports, logs and result folders on F:. Every run creates a new subfolder; old results are preserved.

Input files and reads are referenced directly without copying or uploading. The page binds only to `127.0.0.1`, uses a private per-launch URL, validates Host and request Origin, and loads no external fonts, analytics or web assets. The setup button creates a separate managed Linux environment, downloads core tools and the MOB-suite and PlasmidFinder databases, and checks a PlasmidFinder positive control. On Windows, installation requires the selected WSL distribution to be stored on F:. Environments, databases and caches stay inside that distribution. No administrator password or security-policy change is needed. Existing environments are preserved. Setup logs can be expanded; failed setup can be retried.

## Existing Linux tools

For reconstruction and other Linux tools, expand the environment line in **Run settings**:

- Change **Run analysis using** to **WSL / Linux tools**.
- Set the distribution, usually `Ubuntu-24.04`.
- Click **Find installed environments**. The app suggests the environment with the most matching installed executables. You can choose another from the environment field.
- Click **Check selected tools**. Only tools needed for your current choices are required.

Manual mode remains available for existing installations. Running without marker checks requires explicitly allowing a limited exploratory run. Reports, Excel and the page flag missing evidence; zero candidates does not prove plasmid absence.

Existing selected settings are remembered after a run starts. `app_settings.json` is stored next to the app; remove it to reset settings.

## Workflows

| Workflow | Inputs | What runs |
| --- | --- | --- |
| Screen assembly | Assembly FASTA; optional evidence | Existing candidate screening |
| Reconstruct & check | Assembly FASTA and separate forward/reverse reads | Screening, Unicycler, candidate tracking and eligible junction checks |
| Reuse a completed assembly | Assembly, reads and completed Unicycler folder | Tracking and junction checks using existing assembly/graph |
| Recover a previous run | Original assembly/reads and previous PlasmidTrace workflow folder | Strict recovery of completed saved graphs; not an exact SPAdes resume |
| Investigate candidate genes | Candidate-only FASTA and one or more gene modules | Existing gene adapters; no plasmid classification |

Additional workflows are in **Other workflows**. Paired-read fields appear only when needed. Default resources remain 2 threads and 7 GB. Reconstruction can take hours; the activity indicator does not claim a completion percentage.

## Optional evidence

**Add supporting checks** offers barrnap, MOB-Recon and PlasmidFinder. Automatic mode prepares and enables all three core checks. Manual mode requires installed tools and databases. PlasmidFinder requires its database folder. Expand **Use existing evidence or annotation** to choose existing reports, a matching graph, and GenBank OR GFF3 annotation. Select either a new check or its existing report, not both. Advanced executable paths remain collapsed in **Advanced tool locations**.

Choose Windows-drive files with Browse; paths are converted automatically for WSL. For Linux-only executable or database paths, enter the Linux path or command name directly. Outputs must be in a folder accessible to the Windows launcher, such as F:.

## Candidate gene analysis

This section appears for reconstruction and candidate gene workflows, and stays collapsed until needed. Each module offers **Leave off**, **Use existing results**, or **Run installed tool**; only relevant fields are shown.

- **AMRFinderPlus:** select a TSV or an executable and optional database. `--plus` is retained; stress and virulence elements are not automatically AMR genes.
- **VirulenceFinder:** select a TSV or its Python, BLAST and database folder. The existing example settings remain 0.90 identity and 0.60 coverage.
- **Bioremediation Gene Miner:** provide matching candidate GenBank annotation. For a blank imported TSV, supply the workbook for empty-result verification. Running the module also requires its installed Python, script, DIAMOND reference database, reference metadata and annotation rules.

Reconstruction passes only its selected candidate sequences to enabled modules. Gene-analysis-only mode uses the candidate FASTA you choose. Gene Miner annotation must match candidate DNA exactly or by reverse complement; newly reconstructed circles are not annotated automatically. Screening-only mode leaves cargo modules off even if settings from an earlier workflow remain saved.

## Progress, results and cancellation

The progress panel shows the current state and keeps detailed logs collapsed unless a run fails. Buttons open the report, download Excel, download all results as a ZIP, or open the output folder. ZIP preparation can take time for large runs; original nested outputs and provenance are retained. Excel has standard cell-text limits; original TSV files are authoritative.

**Cancel run** stops only this run and its descendants, including its own Linux process group. It does not shut down WSL or stop unrelated analyses. Partial outputs and logs remain available.

Closing the browser tab does not stop a run. Reopen the private launcher link to return to progress. Stopping the launcher with Ctrl+C cancels the active run. When restarting the launcher, previous results remain in their folders but are not reloaded into a new session.

## Preview status

The existing analysis engine, thresholds, release tag and DOI are unchanged. This app preview has not been published as a new software release. Automated checks exercise the engine, runner, local web API and downloads. Actual Windows file selection, WSL installations and full reconstruction still need target-PC testing. Biological benchmarking remains pending; candidate calls, circularity and read support do not confirm plasmid identity, and gene predictions do not establish activity or plasmid carriage.

## Automatic setup scope

Core setup includes barrnap, MOB-suite, PlasmidFinder, Unicycler, SPAdes and supporting alignment programs. Optional AMRFinderPlus, VirulenceFinder and Gene Miner modules still use their existing advanced configuration; they are not installed by this core setup. Runtime package and database provenance is retained. Setup verifies package compatibility and a PlasmidFinder marker control; it does not benchmark biological performance.
