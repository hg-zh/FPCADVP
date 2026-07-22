# Paper Code Repository

This repository contains the code used to reproduce the simulation studies and real-data analyses in the paper.

Most scripts are written in R. The `RealData/HCPS1200` pipeline also includes Python scripts for data inspection and preprocessing.

## Repository structure

The repository is organized into two main folders:

- `Simulation/`: code for the simulation studies and figures from the simulation section.
- `RealData/`: code for the real-data applications in Section 6 of the paper.

In general, each script should be run from its corresponding subfolder so that relative paths work correctly.

## Simulation

The `Simulation` folder contains three subfolders.

### `Simulation/Demo_Figure2`

This folder provides a quick demo of the proposed method.

- Run `fig2_plots.R` to reproduce the middle and right panels of Figure 2 in the paper.

### `Simulation/Rates_Figure2`

This folder contains the code for the rate/phase-transition experiment.

1. Run `run_phase_transition.R` to generate the results under the `outputs/` folder.
2. Run `plot_phase_transition.R` to reproduce the left panel of Figure 2.

### `Simulation/parallel_run`

This folder contains the Slurm scripts for running the proposed method and the baseline methods on a server.

- Use the `.sbatch` files together with the corresponding parameter settings to launch the jobs.
- After the runs finish, the results are saved in the `Res/` folder.
- Run `analyze_boxplot.R` to generate Figure 3 in the paper.

## Real data

The `RealData` folder contains three subfolders.

### `RealData/Mortality`

This folder contains the code used to reproduce the results in Section 6.1.

- The `Data/` folder contains the preprocessed data.
- Run `01_mean_eigenfunctions_female.R` to reproduce the mean and eigenfunction plots, corresponding to the left panel of Figure 4.
- Run `02_country_sign_heatmap_female.R` to reproduce the right panel of Figure 4.
- Run `03_trajectory_html_female.R` to reproduce the results in Figure 5.

### `RealData/Glucose`

This folder contains the code for the glucose data analysis.

1. Run `run_section62_glucose.R`.
2. Run `plot_section62_glucose.R`.

These two scripts reproduce the results in Figure 6.

### `RealData/HCPS1200`

This folder contains the code for the HCP S1200 application.

1. Download the required data first.
2. Run `01_inspect_data.py`.
3. Run `02_extract_roi_distributions.py`.
4. Run `03_build_transports.R` to construct the preprocessed data objects.
5. Run `run_hcp_eigentransport_fit.R`.
6. Run `make_hcp_eigentransport_modes.R`.

These steps reproduce Figure 7 in the paper.

## Notes

- The scripts write outputs into the corresponding `output/`, `outputs/`, `Res/`, or figure folders inside each subdirectory.
- Some analyses depend on preprocessed data that are already included in the repository, while the `HCPS1200` analysis requires downloading the original data separately before running the pipeline.
