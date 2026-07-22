# Kernel and bandwidth sensitivity analysis

This directory reproduces the Setting 1 sensitivity study for the proposed
transport-FPCA method. It compares Gaussian, Epanechnikov, and quartic temporal
smoothing kernels over joint bandwidth multipliers
\(c_h\in\{0.50,0.75,1.00,1.25,1.50\}\).

## Design

The full study uses 200 Monte Carlo replications at each of 20 cells formed by
\(n\in\{20,50,100,200,400\}\), \(N\in\{3,5,10,20\}\), and \(m=200\).
Within a cell and replication, one simulated data set is shared by all 15
kernel--bandwidth configurations.

The untouched fdapace default is the Gaussian \(c_h=1\) reference. The explicit
Gaussian reference bandwidths are multiplied by \(c_h\) and by a
second-moment normalization factor: 1 for Gaussian, \(\sqrt{5}\) for
Epanechnikov, and \(\sqrt{7}\) for quartic. Thus, each compactly supported
kernel is compared at the same second-moment bandwidth scale as the Gaussian
kernel.

## Placement and requirements

Place this directory directly under the project Code directory:

~~~text
Code/
├── Simulation/
│   └── parallel_run/
│       └── transport_fpca_modular.R
└── sensitivity_final/
    ├── README.md
    ├── run_sensitivity.R
    ├── run_sensitivity.sbatch
    ├── aggregate_results.R
    └── generate_table.R
~~~

The simulation depends on
../Simulation/parallel_run/transport_fpca_modular.R. The driver verifies that
this file has MD5 checksum 4beac3d60a2b3e8850037ea7993abc05. Its corresponding
SHA-256 checksum is
8018f71201c95e8f5321a9faefd65f57f6331162f3a76f6f9bacb24f963c1db0.

Required software is R 4.5.0, fdapace 0.6.0, Bash, and Slurm. The batch script
also uses sha256sum, awk, grep, and squeue. The requested full-run resources
are one node, one task, 50 CPUs, 570 GB of memory, and 120 hours. The partition
and account are intentionally left to the cluster defaults; add the appropriate
sbatch options if the cluster requires them.

From Code/sensitivity_final, run this lightweight preflight before submission:

~~~bash
module load r/4.5.0
Rscript run_sensitivity.R mode=full check_only=TRUE
~~~

## Slurm execution

Create the output directories and submit the full study from
Code/sensitivity_final:

~~~bash
mkdir -p logs results
sbatch \
  --output="$PWD/logs/%x_%j.out" \
  --error="$PWD/logs/%x_%j.err" \
  --export=ALL,R_MODULE=r/4.5.0,SENSITIVITY_CODE_DIR="$PWD",SENSITIVITY_OUTPUT_ROOT="$PWD",SENSITIVITY_LOG_ROOT="$PWD/logs",RUN_TYPE=full \
  run_sensitivity.sbatch
~~~

If the R module has a different name, change R_MODULE. If R 4.5.0 is already
available without an environment module, omit R_MODULE from --export.

After submission, substitute the reported identifier below:

~~~bash
JOB_ID=REPORTED_IDENTIFIER
squeue -j "$JOB_ID"
tail -f "logs/sensitivity_$JOB_ID.out"
~~~

The simulation workload is guarded against execution outside Slurm. A
RUN_TYPE=smoke submission uses two replications at \(n=20\), \(N=3\), and
\(m=200\), but the full supplement table is produced only from RUN_TYPE=full.
When using smoke mode, override --cpus-per-task, --mem, and --time at submission
with smaller values appropriate for the cluster; otherwise the full-run Slurm
directives remain in effect.

## Aggregation and table reproduction

The batch script automatically validates and aggregates a completed run. To
repeat aggregation into a new directory without changing the run output:

~~~bash
RUN_DIR="results/full_REPORTED_IDENTIFIER"
Rscript aggregate_results.R \
  input_dir="$RUN_DIR" \
  outdir="$RUN_DIR/aggregated_manual" \
  mode=full \
  strict=TRUE \
  write_table=FALSE
~~~

Aggregation refuses a nonempty output directory. Regenerate the final LaTeX
table from the four summary CSV files:

~~~bash
Rscript generate_table.R \
  summary_dir="$RUN_DIR/aggregated_manual" \
  table_out="$RUN_DIR/sensitivity_table.tex"
~~~

Table generation requires summary_by_cell.csv,
paired_summary_by_cell.csv, summary_across_cells.csv, and
failure_summary_by_config.csv. It validates the complete design and refuses to
overwrite an existing table.

## Reproducibility and outputs

The random seed for replication \(r\) and cell \((n,N,m)\) is

\[
2000 + 10^6 r + 1000n + 10N + m.
\]

All 15 configurations reuse that cell's simulated data, so paired comparisons
do not introduce independent Monte Carlo variation. The driver records the R
and fdapace versions, dependency checksum, configuration, session information,
and run manifest. Replication files are written atomically.

A full run is stored under results/full_JOB_ID and contains:

- replications/rep_001.rds through replications/rep_200.rds;
- raw_results.rds and raw_results.csv;
- run_config.rds, run_config.csv, run metadata, and sessionInfo.txt;
- batch_manifest.txt and completion markers;
- aggregated validation reports, numerical-status files, cell summaries,
  paired summaries, failure summaries, and sensitivity_table.tex.

Strict aggregation distinguishes structural validity from numerical fit
failures. Compact-kernel failures remain visible in the failure summaries and
table; they do not silently disappear from the planned-fit counts.
