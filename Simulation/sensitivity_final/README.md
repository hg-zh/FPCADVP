# Kernel and bandwidth sensitivity analysis

This directory reproduces the Setting 1 sensitivity study for the proposed
transport-FPCA method. It compares Gaussian, Epanechnikov, and quartic temporal
smoothing kernels over joint bandwidth multipliers
\(c_h\in\{0.50,0.75,1.00,1.25,1.50\}\).

Three procedure settings are supported:

- `latent=FALSE, methodXi=IN`: FPCA and numerical-integration scores on U,
  associated with the dense construction in (13). This is presented first in
  Table S.1, using the representative cell n=100, N=20, m=200. fdapace uses
  trapezoidal integration of centered observations against estimated
  eigenfunctions, a numerical alternative to the uniform-design empirical
  average displayed in the paper.
- `latent=TRUE, methodXi=CE`: FPCA and conditional-expectation scores on
  Z=g^{-1}(U), followed by reconstruction through g, associated with (14).
  This is presented second in Table S.2, using n=100, N=5, m=200.
- `latent=FALSE, methodXi=CE`: the historical direct-U conditional-expectation
  analysis, available as an optional diagnostic. It shares the
  score-prediction principle of (14) but omits its latent-Z transformation;
  it is no longer a displayed sensitivity table in the supplement.

Procedures are run and summarized separately using the same simulation seeds,
and both displayed tables include summaries across all 20 cells. The code uses
the generator's common target baseline norm; subsequent transport projection
can alter the realized subject norms, which are recorded in the diagnostics.
The numerical implementations retain `u_cap=0.95` and `Z_cap=8` (the latter is
active only for the latent procedure). These are regularized implementations
of the constructions discussed in the paper.

## Design

The full study uses 200 Monte Carlo replications at each of 20 cells formed by
\(n\in\{20,50,100,200,400\}\), \(N\in\{3,5,10,20\}\), and \(m=200\).
Within a cell and replication, one simulated data set is shared by all 15
kernel--bandwidth configurations.

Within each procedure, the untouched fdapace kernel/bandwidth default is
the Gaussian \(c_h=1\) reference. The explicit
Gaussian reference bandwidths are multiplied by \(c_h\) and by a
second-moment normalization factor: 1 for Gaussian, \(\sqrt{5}\) for
Epanechnikov, and \(\sqrt{7}\) for quartic. Thus, each compactly supported
kernel is compared at the same second-moment bandwidth scale as the Gaussian
kernel.

## Requirements and paths

Keep the distributed `Simulation/sensitivity_final` and `Simulation/parallel_run`
directories beside each other; no directories need moving. Commands below run
from the public repository root containing `Simulation/` and `RealData/`. If using
the manuscript workspace, first run `cd Code` once. The scripts resolve their
dependency relative to their own location.

The exact simulation environment is R 4.5.0 and fdapace 0.6.0. Preflight checks
both versions and the unmodified transport-FPCA module checksum. The Slurm batch
workflow requires Bash, `sha256sum`, `awk`, `grep`, and `squeue`. The simulation
workload requires a Slurm allocation; preflight and table generation run locally.
The full-run request is one node, one task, 50 CPUs, 570 GB, and 120 hours.
Account, partition, node, installation directories, and module names are not
hardcoded. Supply cluster-specific account/partition options to `sbatch` if your
site requires them. Set `R_MODULE` only if R must be loaded through a module;
otherwise the launcher uses `Rscript` from `PATH`.

```bash
Rscript Simulation/sensitivity_final/run_sensitivity.R mode=full latent=FALSE methodXi=IN check_only=TRUE
Rscript Simulation/sensitivity_final/run_sensitivity.R mode=full latent=TRUE methodXi=CE check_only=TRUE
```

## Full simulation and aggregation

The source-only distribution contains no data files or saved simulation results.
The simulated inputs are generated internally from the seeds and settings below.
Submit the two studies separately:

```bash
sbatch --export=ALL,RUN_TYPE=full,LATENT=FALSE,METHOD_XI=IN Simulation/sensitivity_final/run_sensitivity.sbatch
sbatch --export=ALL,RUN_TYPE=full,LATENT=TRUE,METHOD_XI=CE Simulation/sensitivity_final/run_sensitivity.sbatch
```

The batch script automatically validates and aggregates each completed run.
Default output directories are `results/estimator13_IN` and
`results/estimator14_CE` beneath this directory. They do not depend on a job ID.
The default procedure is estimator (13); specify `LATENT=TRUE,METHOD_XI=CE` for
estimator (14). The optional direct-U CE diagnostic uses
`LATENT=FALSE,METHOD_XI=CE` and writes `results/direct_U_CE`.

An existing output directory is never overwritten. To repeat a completed run,
use a new output root; for example:

```bash
sbatch --export=ALL,SENSITIVITY_OUTPUT_ROOT="$PWD/recomputed_sensitivity",RUN_TYPE=full,LATENT=FALSE,METHOD_XI=IN Simulation/sensitivity_final/run_sensitivity.sbatch
```

For an installation check, replace `RUN_TYPE=full` with `RUN_TYPE=smoke` and
request smaller resources, e.g. `--cpus-per-task=2 --mem=8G --time=01:00:00`.
Smoke mode uses two replications at n=20, N=3, m=200 and all 15 configurations.
Smoke results are written to `results/smoke_estimator13_IN` or
`results/smoke_estimator14_CE`; these do not reproduce the supplement tables.

Repeat strict aggregation into a new output directory, if desired:

```bash
Rscript Simulation/sensitivity_final/aggregate_results.R input_dir=Simulation/sensitivity_final/results/estimator13_IN outdir=Simulation/sensitivity_final/results/estimator13_IN/aggregated_manual mode=full latent=FALSE methodXi=IN strict=TRUE write_table=FALSE
```

Use `estimator14_CE`, `latent=TRUE`, and `methodXi=CE` for the other study.
Aggregation requires the full replication files, configuration, metadata, and
completion marker. The following table-rendering commands use the aggregate
summaries generated by the completed full workflow above.

## Reproduce Tables S.1 and S.2 from the generated summaries

After both full runs and their automatic aggregation have completed:

```bash
mkdir -p Simulation/sensitivity_final/reproduced_tables
Rscript Simulation/sensitivity_final/generate_table.R summary_dir=Simulation/sensitivity_final/results/estimator13_IN/aggregated table_out=Simulation/sensitivity_final/reproduced_tables/sensitivity_table_in.tex
Rscript Simulation/sensitivity_final/generate_table.R summary_dir=Simulation/sensitivity_final/results/estimator14_CE/aggregated table_out=Simulation/sensitivity_final/reproduced_tables/sensitivity_table_latent_ce.tex
```

The renderer validates procedure labels and the full design using
`summary_by_cell.csv`, `paired_summary_by_cell.csv`,
`summary_across_cells.csv`, and `failure_summary_by_config.csv`. It refuses to
overwrite an existing table; use a new `table_out` for a second rendering.

## Randomness and recorded outputs

The cell seed is `2000 + 1000000*r + 1000*n + 10*N + m`. All 15 configurations
reuse that cell's simulated data. Full runs contain 200 replications over 20
cells, hence 60,000 planned fits per procedure, with Gaussian/Epanechnikov/quartic
kernels and bandwidth multipliers 0.50, 0.75, 1.00, 1.25, 1.50. Setting 1 uses
Beta bounds (0.6,3) for both shape parameters, `decayR=1`, random observation
times, 128 transport-grid points, 51 time-evaluation points, and 50 generating
basis functions. Numerical settings are `eps_u=0.05`, `u_cap=0.95`, `Z_cap=8`.

Full runs record replication-level results, run configuration, package versions,
module checksum, fit failures, and validation reports. Strict aggregation
distinguishes structural completeness from numerical fit failures. Compact-kernel
failures remain in the planned-fit counts and table summaries; no failed cell is
silently omitted. Historical execution records are preserved outside the active
public workflow and are not required by any reproduction command.
