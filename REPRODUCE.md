# Reproduction guide

## Layout and software

Run these commands from the **public repository root**, containing `Simulation/`, `RealData/`, and this file. In the manuscript workspace this is `Code/` (`cd Code` once). Do not rename directories or move source files. R entrypoints locate companion code and default inputs relative to their saved location. Quote paths containing spaces. Use explicit relative overrides from the working directory shown here; Rates scripts interpret `outdir` relative to `Simulation/Rates_Figure2/`.

No real-data files, saved application fits, or simulation result datasets are distributed. Simulation scripts generate their own observations. Application data must be obtained separately from the cited sources and supplied in the documented formats. Published graphics in `manuscript_assets/` are reference outputs, not input datasets.

R 4.5.0 and `fdapace` 0.6.0 are the reference simulation environment; the rate and application runs also used R 4.5.1. Fresh sensitivity simulations enforce exactly R 4.5.0 and the original scientific module checksum. Other reference versions are in `dependencies.tsv`. HCP preprocessing needs Python 3, using only its standard library. Optional HTML plots need `plotly`/`htmlwidgets`. Install dependencies explicitly before running scripts, for example in R:

```r
install.packages(c("fdapace", "ggplot2", "frechet", "dplyr", "readr", "tidyr",
                   "plotly", "htmlwidgets"), repos="https://cloud.r-project.org")
# To request the reference fdapace version explicitly:
install.packages("remotes", repos="https://cloud.r-project.org")
remotes::install_version("fdapace", version="0.6.0",
                         repos="https://cloud.r-project.org", upgrade="never")
```

Source installation may need a compiler and the package's dependencies. Then:

```sh
Rscript check_environment.R full=TRUE
```

This check reports versions and fails for missing packages or a different `fdapace`; it never installs or changes software. Other reference-version differences are reported rather than silently corrected.

## Figure and table index

Numbers refer to the revised manuscript/supplement. Output paths below are relative to the public root. Simulation/application result directories are created by the commands; they are not shipped.

| Item | Workflow and procedure | Output |
|---|---|---|
| Figure 1 | Conceptual TikZ diagram; no data | `manuscript_assets/Figure1.tex` |
| Figure 2 left | Rates; direct-U IN, estimator (13); generated Setting 2 data | `Simulation/Rates_Figure2/outputs/estimator13_IN/figures/rate_decay_IN.pdf` |
| Figure 2 middle/right | Demo; direct-U CE, seed 2001, true observed transports | `Simulation/Demo_Figure2/outputs/Fig2_middle_modular.pdf`, `Fig2_right_modular.html`, and `Fig2_right_static.pdf` in the same directory |
| Figure 3 | Comparison; latent CE, estimator (14), versus LQD/Tangent IN; Setting 1 | `Simulation/parallel_run/results/paper_comparison/figures/setting1/Fig_boxplot_byN_methods_m200.pdf` |
| Table 1 | Same comparison, both settings; Monte Carlo mean IMSE | `Simulation/parallel_run/results/paper_comparison/summary/table1_rows.tex` and `table1_display.csv` |
| Figure 4 | Mortality male direct-U IN | `RealData/Mortality/output/estimator13_IN/{mean_IN_M,countries_M}.pdf` |
| Figure 5 | Mortality female Belarus and Spain, direct-U IN | `RealData/Mortality/output/estimator13_IN/{blr_IN_F,esp_IN_F}.pdf` |
| Figure 6 | Glucose latent CE, estimator (14); seed 20 | `RealData/Glucose/output/section62_glucose_seed_20/ind_gluco_seed_20.pdf` |
| Table S.1 | Setting 1 sensitivity, direct-U IN; representative cell: n = 100, N = 20, m = 200 | `Simulation/sensitivity_final/results/estimator13_IN/aggregated/sensitivity_table_in.tex` |
| Table S.2 | Setting 1 sensitivity, latent CE; representative cell: n = 100, N = 5, m = 200 | `Simulation/sensitivity_final/results/estimator14_CE/aggregated/sensitivity_table_latent_ce.tex` |
| Figure S.1 | HCP latent CE; downloaded motor/rest cohort | `RealData/HCPS1200/outputs/figures/hcp_eigentransport_modes.pdf` |
| Table S.3 | Glucose; all 50 eligible subjects / 100 held-out days | `RealData/Glucose/output/section62_glucose_seed_20/holdout_aggregate.csv` |
| Figure S.2 | Mortality female direct-U IN | `RealData/Mortality/output/estimator13_IN/{mean_IN_F,countries_F}.pdf` |
| Figure S.3 | Archived glucose mean/eigenfunction graphic | `manuscript_assets/FigureS3_published.pdf` |

Table 1's generated rows go inside the manuscript's existing table; sensitivity files contain complete table blocks. Table S.3 uses the mean dW and mean squared dW columns of `holdout_aggregate.csv`. The supplied reference PDFs do not constitute a fresh numerical reproduction. Figure S.3's original source fit/exporter is not identified, so no unrelated fit is presented as its exact recreation.

## Figure 1 and Figure 2

With LaTeX and TikZ installed, Figure 1 compiles without any data:

```sh
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=manuscript_assets manuscript_assets/Figure1.tex
```

The full rate workflow is:

```sh
Rscript Simulation/Rates_Figure2/run_phase_transition.R
Rscript Simulation/Rates_Figure2/validate_phase_transition.R
Rscript Simulation/Rates_Figure2/plot_phase_transition.R
```

Settings: 200 replications; n = 20, 30, 40, 50, 70, 100, 150; N = 3, 4, 5, 7, 10, 20; m = 200; seed base 20260323; random observation times. The Setting 2 Beta shape bounds are (3, 4) and (1, 2). The generator uses 50 cosine components, 128 probability points, 51 evaluation times, and `decayR=1`. Fitting uses `latent=FALSE`, `methodXi=IN`, `eps_u=0.05`, `u_cap=0.95`, and `Z_cap=8` (inactive for direct-U fitting). The default display shows n <= 100 and N <= 10; `display=full` retains all 42 cells. The reference exponent is -(a-1)/(3a+2), with a = 2. Archived CE results are optional, not a prerequisite. To preserve a previous run, use `outdir=outputs/new_IN` consistently for all three commands; `cores=5` controls computation in the driver.

The demonstration is independent of the rate experiment:

```sh
Rscript Simulation/Demo_Figure2/fig2_plots.R
```

The demonstration uses n = 100, N = 4, m = 200, 50 cosine components, seed 2001, subject 9, the Setting 2 Beta bounds, and direct-U CE reconstruction of true observed transports. Default outputs are under `Simulation/Demo_Figure2/outputs/`: a middle-panel PDF, interactive right-panel HTML with its dependency directory, and a deterministic static right-panel PDF from the same fit. To select another destination, append `outdir=results/figure2_demo`; this override is relative to the current working directory. The exact published panels are retained as `manuscript_assets/Figure2_*_published.*`. Original visual-export settings were not recorded, so the static companion is not claimed to be a pixel-identical rendering of the archived right panel.

## Figure 3 and Table 1

Run the two site-neutral Slurm jobs from the repository root, then summarize after both finish:

```sh
sbatch Simulation/parallel_run/run_proposed_latent_CE.sbatch
sbatch Simulation/parallel_run/run_benchmarks_IN.sbatch
Rscript Simulation/parallel_run/summarize_squared_imse.R
Rscript Simulation/parallel_run/analyze_boxplot.R setting=1
Rscript Simulation/parallel_run/analyze_boxplot.R setting=2
```

Figure 3 uses Setting 1; the last command produces an additional Setting 2 plot. Each job requests one node, one task, 50 CPUs, 300 GB of memory, and 120 hours, and runs both settings. Scripts use `Rscript` on PATH. Supply your cluster's account/partition options to `sbatch` if required; optionally set `R_MODULE` to its R module name. No account, node, installation location or partition is hardcoded. Submitting from another directory requires `SIMULATION_CODE_DIR`. All full Rscript calls are visible in the launchers and can be run directly on an appropriately resourced non-Slurm machine with `cores=...`.

Launchers refuse to overwrite existing study directories. For a second study, pass `--export=ALL,RESULTS_ROOT="$PWD/recomputed_comparison"` to both `sbatch` calls, then `results_dir=recomputed_comparison` to the summarizer and plotting commands. No editing of source paths is needed.

Both settings use 200 replications, seed base 2000, random observation times, n = 20, 50, 100, 200, 400, N = 3, 5, 10, 20, and m = 200. The generator uses `decayR=1`, 50 cosine components, 128 probability points, and 51 evaluation times. The Beta shape bounds are (0.6, 3) for both parameters in Setting 1, and (3, 4) and (1, 2) in Setting 2. Each cell uses seed `2000 + 1000000*rep + 1000*n + 10*N + m`. Proposed fitting uses `latent=TRUE`, `methodXi=CE`, `eps_u=0.05`, `u_cap=0.95`, and `Z_cap=8`. Benchmarks use `methodXi=IN`, ten cubic B-spline basis functions, and no derivative smoothing.

The summarizer validates all 4,000 planned records per method and setting, together with the saved estimator metadata; failed or missing fits cause failure. It uses the `imse_sq*` fields: the subject average of the time integral of squared Wasserstein distance. At each time, the absolute transport difference is integrated over the common support and then squared; both integrations use trapezoidal quadrature. Table percentages are `100 * (1 - mean proposed IMSE / mean benchmark IMSE)`. Figure 3 and Table 1 display N = 3, 5, 10; full summaries retain N = 20.

## Sensitivity Tables S.1 and S.2

Full simulations require Slurm, R 4.5.0, and fdapace 0.6.0:

```sh
Rscript Simulation/sensitivity_final/run_sensitivity.R mode=full latent=FALSE methodXi=IN check_only=TRUE
Rscript Simulation/sensitivity_final/run_sensitivity.R mode=full latent=TRUE methodXi=CE check_only=TRUE
sbatch --export=ALL,RUN_TYPE=full,LATENT=FALSE,METHOD_XI=IN Simulation/sensitivity_final/run_sensitivity.sbatch
sbatch --export=ALL,RUN_TYPE=full,LATENT=TRUE,METHOD_XI=CE Simulation/sensitivity_final/run_sensitivity.sbatch
```

The launcher validates and aggregates completed runs automatically. See [its README](Simulation/sensitivity_final/README.md) for resource requests, smoke mode, strict reaggregation, and regenerating LaTeX tables from the summaries produced by these commands. No sensitivity data cache is included.

Both sensitivity studies use Setting 1, 200 replications, and the comparison study's 20-cell grid and seeds. Each dataset is shared across 15 configurations: Gaussian, Epanechnikov, and quartic kernels combined with bandwidth multipliers 0.50, 0.75, 1, 1.25, and 1.50, with second-moment normalization. Table S.1 uses `latent=FALSE`, `methodXi=IN`, and the representative cell n = 100, N = 20, m = 200. Table S.2 uses `latent=TRUE`, `methodXi=CE`, and n = 100, N = 5, m = 200. All-cell summaries retain planned-fit counts and failures.

## Application setup

Readers must obtain the data themselves. Data availability does not imply that current downloads exactly match the historical prepared inputs. The local application READMEs specify file formats, ordering, cohort selection and preprocessing requirements. Data and fit objects are excluded from this distribution.

**Mortality (Figures 4, 5, and S.2):** [Human Mortality Database](https://www.mortality.org/). Prepare `Data_Female.RData`, `Data_Male.RData`, `Data_raw_female.Rdata`, `Data_raw_male.Rdata` in `RealData/Mortality/Data/`, or set `MORTALITY_DATA_DIR`. Required objects, dimensions and country ordering are specified in the [Mortality README](RealData/Mortality/README.md).

```sh
Rscript RealData/Mortality/04_rerun_IN_comparison.R
```

Mortality settings: direct-U IN scores, `kappa=500`, three principal components, 33 countries/regions, annual observations from 1983 through 2018, 401 age-grid points over 20--100, `u_cap=0.999999`, no Z cap, and smoothed mean/covariance estimation. There is no random sampling. Figure 4 uses males; Figure 5 shows female Belarus and Spain; Figure S.2 uses females. A completed default output is protected; supply a new positional output-directory argument to rerun. The historical raw-HMD-to-quantile script is not included; exact refitting requires prepared inputs matching the documented schema and original processing.

**Glucose (Figure 6 and Table S.3):** Hall et al. (2018), [PLoS Biology, e2005143](https://doi.org/10.1371/journal.pbio.2005143), and the [functional-diabetes resource](https://prediccion.shinyapps.io/functional_diabetes/) cited in the paper. Supply `RealData/Glucose/Data/data_1.csv` with columns `id`, `time`, and `gl`, or set `GLUCOSE_DATA_FILE`. Preserve row order. See the [Glucose README](RealData/Glucose/README.md).

```sh
Rscript RealData/Glucose/run_section62_glucose.R
Rscript RealData/Glucose/run_section62_glucose_baselines.R
Rscript RealData/Glucose/plot_section62_glucose.R
Rscript RealData/Glucose/summarize_section62_glucose.R
```

Glucose settings: seed 20, proposed latent CE, and CE benchmarks. Two days are held out per eligible subject, defined as having at least six observed days within the first ten-day window. The support is 50--160, with 201 support/probability points and time grid 0, 0.25, ..., 10. Other settings are bandwidth 4, the Epanechnikov kernel, at least 30 readings per day, `spar=0.5`, density threshold 0.01, `kappa_original=3200`, `u_cap=0.95`, and `Z_cap=8`. Timestamp parsing explicitly uses America/Chicago and UTC date conversion, preserving the historical split. Figure 6 selects the two subjects with the smallest proposed held-out error; Table S.3 includes all 50 eligible subjects. Raw-source-to-CSV assembly is not provided. Figure S.3's unidentified original descriptive fit/exporter is distinguished from this predictive workflow; the archived graphic is supplied without claiming exact recomputation.

**HCP (Figure S.1):** [HCP Young Adult/S1200](https://www.humanconnectome.org/study/hcp-young-adult). Readers download the motor/rest ZIP files through the source's access process. Follow the [HCP README](RealData/HCPS1200/README.md) for filenames/layout and external-root options. No data or cohort list is supplied.

```sh
python3 RealData/HCPS1200/01_inspect_data.py
python3 RealData/HCPS1200/02_extract_roi_distributions.py --max-subjects 0 --sample-per-roi 1500 --quantile-grid-size 101 --seed 20260304
Rscript RealData/HCPS1200/03_build_transports.R
Rscript RealData/HCPS1200/run_hcp_eigentransport_fit.R
Rscript RealData/HCPS1200/make_hcp_eigentransport_modes.R
```

Use Python `--root` and R `root=` consistently for external input/output locations. `--max-subjects 0` processes the full downloaded cohort; the default of 20 is only a pilot. The workflow uses four regions, six runs, sampling seed 20260304, up to 1,500 grayordinates per region, 101 quantiles, and 201 transport-grid points. Fitting uses latent CE, maximum-observed-norm normalization, `u_cap=0.95`, `Z_cap=8`, `maxK=3`, and `FVEthreshold=0.999`; the plotted times are 0.40 and 0.925. The manuscript used 209 subjects; arbitrary current downloads, even with the same count, need not reproduce its exact cohort. Inspection/extraction record the actual downloaded cohort and sampled locations. The original cohort manifest is not distributed, and full raw-data HCP reproduction has not been verified in this revision.

## What was verified

[VERIFICATION.md](VERIFICATION.md) records executed checks. Path changes did not trigger full reruns of the long simulations. Historical research outputs are kept outside the public code. Remaining preprocessing and archived-visual limitations are stated above rather than described as completed raw-data reproductions.
