# Automatic browser app preview validation

Based on repository commit `79428741f110751f3d62e69f16864912c45f18dd`. Interface preview 0.4.0 retains engine 1.0.0 classifications and thresholds. The core runner now accepts an explicit MOB-suite database directory.

80 automated tests pass: original engine tests, runner and cancellation tests, local HTTP security and downloads, runtime selection, missing-evidence guards, database readiness and positive-control requirements. Installer orchestration tests simulate external programs and cover success and failed controls; they are not live package-installation tests.

The previous uploaded 112-contig run was reassessed: all marker checks and coverage were unassessed, and the new assessment correctly flags its limitations. Excel, HTML and JSON exports carry assessment warnings.

Python and JavaScript syntax and whitespace checks pass. A Chromium download was truncated, so rendered browser QA is not claimed. Actual Windows file selection, wsl.exe, live dependency/database installation and full reconstruction require target-PC verification. No biological sensitivity or specificity claim is made.

The citation metadata, DOI and published release tag are unchanged. This is a local preview, not a published release.

Windows launcher hotfix: execute a UTF-8 shell script file through WSL instead of passing multiline shell commands through Windows quoting. A regression test checks preservation of shell variables and quotes. Actual wsl.exe verification remains pending.

Simplified interface: 80 Python tests and 3 Node UI handoff checks pass. The latter execute the setup-completion handler for ready, failed and cancelled setup, preserving the selected analysis and preventing launch after unsuccessful setup. Browser rendering and live Windows installation remain unverified.
