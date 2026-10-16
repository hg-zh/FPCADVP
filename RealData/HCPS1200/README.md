# HCP S1200 distributional fMRI application

In the manuscript workspace, run `cd Code` once before these commands. The public repository already starts at this level; do not move or rename its directories.

No downloaded data, cohort list, or processed fits are distributed. Obtain authorized access to the HCP S1200 data resource cited in the paper. Python 3 is required; both Python scripts use only the standard library. R packages are `dplyr`, `readr`, `tidyr`, `fdapace`, and `ggplot2`.

## Input layout

By default, inputs and outputs are relative to this application directory:

```text
RealData/HCPS1200/
  HCP_Data/
    MotorTask/<subject>_3T_tfMRI_MOTOR_preproc.zip
    RestingState/<subject>_3T_rfMRI_REST_fix.zip
```

Keep the original ZIP contents and filenames. The code reads CIFTI `.dtseries.nii` files for motor runs `tfMRI_MOTOR_LR`, `tfMRI_MOTOR_RL` and resting runs `rfMRI_REST1_LR`, `rfMRI_REST1_RL`, `rfMRI_REST2_LR`, `rfMRI_REST2_RL`. Only subjects present in both ZIP collections are used, in sorted identifier order. Inputs may instead be held in an external directory with the same `HCP_Data/` subtree; use Python `--root /path/to/hcp-work` and R `root=/path/to/hcp-work` consistently. Code is always loaded beside the entrypoint, independent of this data/output root.

The manuscript cohort has 209 subjects. Its original exact subject manifest is not present in this distribution; exact reproduction requires that cohort and the original input releases. Using an arbitrary current collection of HCP files will not establish an exact reproduction, even if it has 209 subjects. The inspection command records the supplied cohort, and extraction records the sampled grayordinates and run-level status. No data are downloaded automatically.

## Commands from the public repository root

```sh
python3 RealData/HCPS1200/01_inspect_data.py
python3 RealData/HCPS1200/02_extract_roi_distributions.py --max-subjects 0 --sample-per-roi 1500 --quantile-grid-size 101 --seed 20260304
Rscript RealData/HCPS1200/03_build_transports.R
Rscript RealData/HCPS1200/run_hcp_eigentransport_fit.R
Rscript RealData/HCPS1200/make_hcp_eigentransport_modes.R
```

These entrypoints also work from another directory when given their script paths. The extractor's default is a 20-subject pilot; **`--max-subjects 0` is required to process the full supplied cohort**. Default filenames and statistical settings above are explicit; intermediate directories are created automatically. For external files, pass `--root /path/to/hcp-work` to both Python commands and `root=/path/to/hcp-work` to all three R commands. Relative R path overrides are relative to that root. The R builder accepts `input_quant=...`, `output_quant_norm=...`, `output_transports=...`, `output_mean_quant=...`, `output_objects=...`, `output_norm_stats=...`, `log_file=...`, and `xgrid_len=201`. The fit script accepts `input_transport_rds=...`, `output_fit_rds=...`, and optional `method_script=...`; the plot accepts `fit_rds=...` and `out_pdf=...`.

## Settings and outputs

The fixed sampling seed is 20260304, with at most 1,500 sampled grayordinates per region. The regions are left cortex, right cortex, cerebellum, and subcortical structures, defined using CIFTI brain-model metadata. Each grayordinate contributes its temporal standard deviation, with population divisor (number of frames). The distribution is summarized at 101 equally spaced probabilities and rescaled per region. Four rest times are 0.10, 0.30, 0.50, 0.70; motor times are 0.85 and 1.00. Transport grids contain 201 points. Fitting uses `latent=TRUE`, `methodXi="CE"` (estimator (14)), maximum-observed-norm normalization, `u_cap=0.95`, `Z_cap=8`, `maxK=3`, and `FVEthreshold=0.999`. The figure evaluates times 0.40 and 0.925 and the leading mean-plus/minus-two-standard-deviation modes.

| Stage | Outputs relative to data/output root |
| --- | --- |
| Inspect | `data_processed/hcp_zip_inventory.csv`, `hcp_data_inventory.json` |
| Extract | `data_processed/hcp_roi_quantiles_raw.csv`, `hcp_run_inventory.csv`, `hcp_roi_sampling.json` |
| Build | `data_processed/hcp_distribution_quantiles_normalized.csv`, `hcp_transports_long.csv`, `hcp_mean_quantiles.csv`, `hcp_roi_norm_stats.csv`, `hcp_transport_objects.rds` |
| Fit | `outputs/hcp_transport_fpca_fits.rds` |
| Supplement HCP mode figure | `outputs/figures/hcp_eigentransport_modes.pdf` |

The raw quantile CSV can be supplied as an intermediate input with columns `subject_id, condition, run, t, roi, p, q_raw`. This bypasses extraction but still requires the original cohort and sampling provenance for exact reproduction. To update the manuscript graphic after validation, copy the final PDF to `../paper/Fig/hcp_eigentransport_modes.pdf` when working from the manuscript workspace's `Code/` directory. Each stage records a log under `logs/`. Rerunning overwrites outputs with the same names.

Path resolution, argument parsing, and explicit missing-input failures have been checked on the distributed code. The full raw HCP workflow has not been rerun in this revision because the cohort inputs and original manifest are unavailable.
