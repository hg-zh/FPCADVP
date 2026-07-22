# Mortality Reproduction Scripts

These scripts reproduce the female mortality figures using the read-only implementation in `proposed method new code/transport_fpca_modular.R`.

## Script map

- `00_run_all.R`
  Builds all outputs in one pass.
- `01_mean_eigenfunctions_female.R`
  Replaces `old code/Final_mean_eig.R` and targets `old figure/mean_F.pdf`.
- `02_country_sign_heatmap_female.R`
  Replaces `old code/Country_plot_new.R` and targets `old figure/countries_F.pdf`.
- `03_trajectory_html_female.R`
  Replaces `old code/Individual_illu.R` and targets the legacy HTML files in `old figure/`.

## How to run in RStudio

1. Open the `Mortality` folder in RStudio.
2. Source `00_run_all.R`, or source any individual script.
3. Edit the parameter block at the top of the script before re-running if you want different paths, output names, figure sizes, colors, fonts, legend settings, method settings, or sex.

## Sex-specific filenames

- Set `sex <- "Female"` to save outputs such as `mean_F.pdf`, `countries_F.pdf`, `bel_F.html`, and `blr_F.html`.
- Set `sex <- "Male"` to save outputs such as `mean_M.pdf`, `countries_M.pdf`, `bel_M.html`, and `blr_M.html`.

## Notes on matching the legacy outputs

- The helper file computes the observed mortality transports from the data and then calls `fit_transport_fpca()` from the read-only `proposed method new code` folder whenever the proposed method is needed.
- The mean/eigenfunction plot defaults to the legacy FVE labels from the old plotting script so the legend matches the saved reference figure. Set `use_legacy_fve_labels <- FALSE` in `01_mean_eigenfunctions_female.R` if you want the labels recomputed from the fitted eigenvalues instead.
- The country heatmap script defaults to the row order seen in the saved legacy PDF so the visual layout matches the reference more closely. Remove or edit `ordered_country_codes` in `02_country_sign_heatmap_female.R` if you prefer the rows to be re-sorted directly from the recomputed scores.
- In the original `old code/Individual_illu.R`, `bel_F.html` was saved from `Pos_fig[[3]]`, which is Spain (`ESP`), not Belgium (`BEL`). The new trajectory script mirrors that archived mapping by default with `ESP -> bel_*.html` and `BLR -> blr_*.html`.
- If you prefer true country-code naming instead of the archived filenames, edit `legacy_output_country_map` in `03_trajectory_html_female.R` and `trajectory_file_map` in `00_run_all.R`.
