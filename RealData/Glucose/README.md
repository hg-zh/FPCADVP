# Glucose application

In the manuscript workspace, run `cd Code` once before these commands. The public repository already starts at this level; do not move or rename its directories.

All entrypoints resolve helpers, inputs, and default outputs from their own location, not the working directory. Run the commands below from the public repository root (containing `RealData/` and `Simulation/`), or use the full script path from another directory. R packages required are `fdapace`, `frechet`, and `ggplot2`; `splines` and `grid` come with R. `rstudioapi` is optional.

## Input setup

The source measurements are the Hall et al. (2018) glucose data cited in the paper, also described at the functional-diabetes data resource cited in the Data Availability statement. The prepared CSV and saved application fits are not distributed; readers obtain the data themselves. Supply `RealData/Glucose/Data/data_1.csv`, or set the environment variable `GLUCOSE_DATA_FILE` to its location for the analysis, baseline, and summary commands. Required columns are:

| Column | Format |
| --- | --- |
| `id` | Subject identifier, consistent across observations |
| `time` | Timestamp text, `YYYY-MM-DD HH:MM:SS` |
| `gl` | Glucose reading in mg/dL; nonnumeric values become missing |

The original file contains 96,412 rows and 52 subjects. Preserve its row order: subject order is the first-occurrence order, and the first recorded date defines each subject's day origin. An optional `X` index column is accepted. The loader removes incomplete rows, matching the original analysis. Timestamp parsing is explicitly `America/Chicago`, with conversion to calendar dates in `UTC`, which reproduces the historical computation. The optional `GLUCOSE_TIMEZONE` override is for other datasets; changing it changes the analysis. Raw-source-to-CSV assembly is not supplied, so exact reproduction requires this prepared CSV.

## Commands

```sh
Rscript RealData/Glucose/run_section62_glucose.R
Rscript RealData/Glucose/run_section62_glucose_baselines.R
Rscript RealData/Glucose/plot_section62_glucose.R
Rscript RealData/Glucose/summarize_section62_glucose.R
```

For external inputs, set `GLUCOSE_DATA_FILE=/path/to/data_1.csv` before those commands; no changes to code or directory layout are needed. Output files are in `RealData/Glucose/output/section62_glucose_seed_20/`. The plotting script reconstructs its output path locally instead of using absolute paths embedded in older saved objects. Rerunning these commands overwrites files in that output directory.

## Settings and manuscript mapping

The seed is 20. Subjects with at least six study days have two days held out, without replacement; all available days within the first ten-day window are used. The published split has 50 eligible subjects and 100 held-out days. Proposed fitting uses `latent=TRUE`, `methodXi="CE"` (estimator (14)), `kappa_original=3200`, `u_cap=0.95`, `Z_cap=8`, and `error=FALSE`. The original glucose support is 50--160 with 201 grid points; the probability grid also has 201 points; fitting times are 0--10 in steps of 0.25. Fréchet regression uses bandwidth 4 and the Epanechnikov kernel; density threshold 0.01, at least 30 points per day, and spline smoothing parameter 0.5. The benchmark scripts use CE scores, ten cubic B-spline basis functions, and the same training split and mean transports. See the saved `config` objects for the complete executed settings.

| Manuscript item | Output and producing command |
| --- | --- |
| Main Figure 6 | `ind_gluco_seed_20.pdf`, final plot from `plot_section62_glucose.R` |
| Supplement aggregate glucose prediction table | `holdout_aggregate.csv`, from `summarize_section62_glucose.R` |
| Full prediction audit | `holdout_per_day.csv`, `holdout_per_subject.csv`, `holdout_aggregation_audit.txt`, from summary command |
| Saved proposed and benchmark fits | `section62_results_seed_20.rds`, `baseline_results_seed_20.rds`, first two commands |
| Supplement glucose mean/eigenfunction figure (`Mean_gl`) | Archived manuscript graphic; its original plotting script and exact source fit are not identified in the distributed workflow |

Figure 6 selects the two subjects with the smallest proposed-method mean held-out error for illustration. The aggregate table includes every eligible subject. The summary computes the square of the integrated absolute transport difference (integrating first, then squaring), and also reports the unsquared criterion. Stored means and predictions can be plotted without the raw CSV; eligibility and split validation in the summary still require the CSV. The archived `Mean_gl` visual is not silently regenerated from the seed-20 predictive fit, because equivalence to its original full-data fit has not been established.
