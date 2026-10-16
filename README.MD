# Functional Principal Component Analysis for Distribution-Valued Processes

The public repository root contains `Simulation/` and `RealData/`. Keep this layout as distributed. In the manuscript workspace, enter `Code/` once before running the commands below. No source-path edits are needed.

Start with [REPRODUCE.md](REPRODUCE.md): it maps every manuscript figure and table to its commands, inputs, estimator settings, and outputs. Reference dependency versions are in [dependencies.tsv](dependencies.tsv).

```sh
Rscript check_environment.R full=TRUE
Rscript Simulation/Demo_Figure2/fig2_plots.R
```

The second command runs a self-contained simulation demonstration. Full studies and application commands are documented in the guide. Real-data inputs, saved fits, and simulation result datasets are not bundled. Readers supply application data from the cited sources; scripts generate simulation data themselves. Application READMEs describe required filenames, formats, and external input locations. Missing inputs produce setup messages instead of requiring code edits.

The portable distribution contains active source, documentation, and published reference graphics. Historical execution logs, generated result archives, and empirical-data fit objects are excluded. See [VERIFICATION.md](VERIFICATION.md) for checks actually performed. The guide distinguishes verified code execution from unavailable historical preprocessing or visual-export provenance.
