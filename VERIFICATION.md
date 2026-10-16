# Reproduction verification, 15 October 2026

The public distribution uses the existing root layout: `Simulation/`,
`RealData/`, and the root README. It contains no input datasets, simulation
result datasets, or saved application fits. Static published graphics are
identified as reference outputs.

## Executed checks

- All 17 active simulation R scripts and 19 real-data R scripts parsed; all
  six Slurm launchers passed Bash syntax checks. The scientific transport/FPCA
  modules were unchanged by this portability revision.
- The dependency checker passed on R 4.5.1 with fdapace 0.6.0 and the listed
  reference package versions, both from the code root and an unrelated working
  directory. `source()` use and paths containing spaces were also checked.
- A relocated source-only real-data tree, with spaces in its directory name,
  was exercised from another working directory. All 12 R entrypoints reached
  the appropriate missing-input/setup messages; both HCP Python command-line
  parsers worked. These were setup checks, not full refits without data.
- A paired Setting 1 simulation smoke run (one replication, n=20, N=3, m=200)
  completed for proposed latent CE and both IN benchmarks with finite squared
  IMSE. The full comparison settings are preserved in the batch scripts.
- Reprocessing the separately held research results validated all 4,000
  planned records per method and setting, regenerated Table 1, and plotted
  both settings. All 8,400 Figure 2 IN records passed standalone validation.
  Both sensitivity tables regenerated from the existing aggregate summaries.
  Those result datasets are not included in the public code.
- The Figure 2 demonstration generated its middle PDF, interactive right
  panel, and static companion from a relocated invocation. The original
  published right-panel export remains a separately identified reference.
- Glucose plotting regenerated Figure 6 from separately held saved fits in a
  relocated tree, demonstrating that old embedded output paths are ignored.
  Explicit timezone parsing preserved the recorded subject list, study days,
  and seed-20 holdout split exactly before the input CSV was removed.

The Figure 2 cache validation retains the original driver checksum and
recognizes the exact old/new pair differing only by portable filename decoding.
The scientific-module checksum still must match. New runs record current source
checksums. Original records were not relabeled as fresh simulation runs.

## Scope and remaining setup

The full 200-replication studies were not rerun solely to test path edits.
Fresh sensitivity simulations require R 4.5.0 exactly; the local R 4.5.1
preflight correctly rejects that different version. Slurm launchers were
syntax-checked, not submitted to an unspecified new cluster.

Application data must be downloaded/prepared by readers. No full raw-data
application refit was claimed after removal of those inputs. The historical
raw-HMD-to-quantile and glucose-source-to-CSV scripts, exact HCP cohort manifest,
and original descriptive glucose figure exporter are not supplied. The guide
states the required formats and identifies these limits. A new downloaded HCP
cohort may legitimately differ from the manuscript's historical cohort.

## Release checks

`tools/build_release.py` uses an explicit source/document/static-asset allowlist,
rejects data files and personal installation paths in distributed text, and
includes `SHA256SUMS` covering every other ZIP member. Its layout, deterministic
archive bytes, checksum manifest, rejection behavior, and overwrite protection
were checked. The actual archive was extracted into a fresh directory containing
spaces: all 70 file hashes passed, no datasets were present, all 37 R files and
three Python files parsed, and all six Slurm scripts passed syntax checks.
Dependency checks worked from the archive root and an unrelated directory;
the self-contained Figure 2 demonstration ran and the Figure 1 TikZ source
compiled there. The release is intended for publication; preparing it locally
does not update the remote repository.
