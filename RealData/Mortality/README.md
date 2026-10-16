# Mortality: direct-U reconstruction with integration scores

In the manuscript workspace, run `cd Code` once before these commands. The public repository already starts at this level; do not move or rename its directories.

The current mortality figures use the numerical implementation associated with estimator (13): `latent=FALSE`, `methodXi="IN"`. The 33 countries/regions have observations at every calendar year from 1983 through 2018. The complete regular grid permits numerical score integration. The local implementation in `transport_fpca_modular.R` is unchanged.

## Reproduce the analysis

From the public repository root (the directory containing `RealData/` and `Simulation/`), run:

```sh
Rscript RealData/Mortality/04_rerun_IN_comparison.R
```

The driver also works from another working directory when given its absolute path. It requires R packages `fdapace` and `ggplot2`; the optional HTML scripts additionally require `plotly` and `htmlwidgets`. The verified run used local R 4.5.1 and fdapace 0.6.0. No cluster is needed.

The input data are excluded from this distribution. Supply the following four prepared files in `RealData/Mortality/Data/`, or set `MORTALITY_DATA_DIR` to the directory containing them:

- `Data_Female.RData` and `Data_Male.RData` contain `Data_q`.
- `Data_raw_female.Rdata` and `Data_raw_male.Rdata` contain `Data_female` and `Data_male`, respectively, including country codes.

`Data_q` must be a list of 33 countries/regions, each containing 36 annual entries ordered from 1983 to 2018. Each entry is `list(x=..., y=...)`, with `x=seq(0,1,length.out=401)` and the corresponding nondecreasing age-at-death quantiles `y` on ages 20--100. The raw objects must be data frames with a `country` column; `unique(country)` in its existing row order must match the outer `Data_q` order. The driver does not use the additional `Data_den` and `Data_p` objects. Country order and original row order must be preserved.

The source observations are Human Mortality Database life tables, as cited in the paper. The historical raw-life-table-to-prepared-quantile script is not present in this distribution. Therefore exact refitting requires the prepared files described above; downloading current life tables alone is not an exact reproduction recipe. No automatic download or replacement smoothing is performed.

For inputs stored outside the repository, no source edits are needed:

```sh
MORTALITY_DATA_DIR=/path/to/prepared-mortality Rscript RealData/Mortality/04_rerun_IN_comparison.R
```

The driver checks required files before creating output. Missing inputs produce an explicit setup message.

The default destination is `output/estimator13_IN/` within this Mortality folder. A completed run is protected against accidental overwrite; to rerun, supply a new destination:

```sh
Rscript RealData/Mortality/04_rerun_IN_comparison.R results/mortality_IN_check
```

## Outputs and manuscript figures

| Output | Content |
| --- | --- |
| `mean_IN_M.pdf`, `countries_M.pdf` | Male direct-U mean/eigenfunctions and empirical transport heatmap (main Figure 4) |
| `blr_IN_F.pdf`, `esp_IN_F.pdf` | Female IN trajectories for Belarus and Spain (main Figure 5) |
| `mean_IN_F.pdf`, `countries_F.pdf` | Female direct-U mean/eigenfunctions and empirical transport heatmap (supplement) |
| `mortality_F_fits.rds`, `mortality_M_fits.rds` | Matched IN and CE fits, reconstructed surfaces, observed transports, and former latent descriptive fit |
| `comparison_by_country.csv`, `comparison_summary.csv` | Country/sex differences and sex-level summaries |
| `fitted_eigenvalues_FVE.csv`, `positive_covariance_spectrum.csv` | Retained eigenvalues, cumulative FVE, and its denominator |
| `CE_IN_covariance_invariance.csv` | Verification that changing the score method leaves the mean, covariance, eigenfunctions and eigenvalues unchanged |
| `fit_diagnostics.csv`, `descriptive_old_latent_vs_direct_U.csv` | Clipping/fit diagnostics and separate comparison against the former latent descriptive plots |
| `empirical_signed_mass_by_country_year.csv`, `displayed_country_mapping.csv` | Observed signs and explicit country-to-file mapping |
| `configuration.*`, `sessionInfo.txt`, `input_code_md5.csv`, `run_complete.txt` | Settings, runtime, input/code provenance and validation |

The original `bel_F.html` was generated from the `ESP` country (Spain), not `BEL` (Belgium); the updated figure therefore uses Spain and the filename `esp_IN_F.pdf`. `blr_F.html` corresponds to `BLR` (Belarus). The illustrative countries are preserved; only the incorrect Belgium label is corrected. Historical outputs are preserved outside the public code distribution. Generated fits and data are not bundled; readers must supply the inputs described above.

## Fixed settings and interpretation

The matched comparison changes only `methodXi` between `CE` and `IN`. Both use direct-U FPCA (`latent=FALSE`), `kappa=500`, three components, smoothed mean/covariance estimation, and the supplied annual grid. The executed dense-design fdapace defaults use the Epanechnikov kernel and mean/covariance bandwidths 1.75/3.5 years. Three components are fixed explicitly for both sexes (the automatic 99% threshold would select only two for males). As in the archived trajectory script, reconstructions omit the estimated mean and retain the baseline upper-clamping/cumulative-maximum adjustment, with no further reconstructed-row projection. IN scores are fdapace's trapezoidal integrals of centered observations; they numerically implement score integration on this fixed grid rather than the equal-weight random-design sum in equation (12).

The inherited `u_cap=0.999999` clips 1 of 1,188 female scalar observations and 97 of 1,188 male observations under `kappa=500`. This behavior is identical in the paired fits and is retained to isolate the score-method change. All numerical summaries use the full 401-point age grid on 20--100 and all 36 years. Surface meshes alone are thinned for compact vector PDFs.

The mean/eigenfunction figures now describe the same rescaled direct-U process as the trajectory fit. The old descriptive plots instead used `latent=TRUE`, `kappa=2000`, and IN scores; that separate change of target/normalization is audited separately, not attributed to the CE-to-IN change. Cumulative FVE labels come from the fitted covariance spectrum (`fpca_fit$cumFVE`), not hard-coded legacy percentages or normalization by only the retained eigenvalues. The male percentages are 96.33%, 99.00%, 99.53%; the female percentages are 93.84%, 98.27%, 99.55%. The empirical sign heatmaps do not use fitted FPCA scores and therefore are unaffected by CE versus IN.

## Matched comparison

Integrals use trapezoidal quadrature in age and calendar time. The mean absolute and root mean squared (RMS) differences divide the respective integrated discrepancies by the age/time domain size before reporting age units (years). Relative RMS differences divide by the observed transport-minus-identity RMS, not by the transport itself; the latter could misleadingly make differences look negligible. These are in-sample comparisons, not prediction validation.

| Sex | Mean absolute IN-CE difference (years) | RMS difference (years) | RMS difference / observed deviation | Maximum pointwise difference (years) | IN / CE in-sample RMSE (years) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Female | 0.0668 | 0.1033 | 4.30% | 0.7481 | 0.6385 / 0.6326 |
| Male | 0.0765 | 0.1124 | 2.91% | 0.5353 | 0.9046 / 0.9147 |

The score change has modest overall effects but is not identical, and IN does not uniformly improve fit. For the displayed female Belarus and Spain examples, RMS differences are 0.2547 and 0.1274 years (5.35% and 4.38% of observed deviations). Mean/covariance/eigenfunction estimates are exactly unchanged by the matched score-method switch. Relative to the old latent descriptive analysis, sign-aligned eigenfunction inner products are 0.99999, 0.99992, 0.99991 for females and 0.99801, 0.98192, 0.95496 for males.

## Optional figure-by-figure scripts

`00_run_all.R` and `01_mean_eigenfunctions_female.R`--`03_trajectory_html_female.R` can be run from any working directory by giving the script path. They now default to direct-U IN fitting, fitted cumulative FVE, and true country-code HTML names (`esp_F.html`, `blr_F.html`). Set `sex <- "Male"` in their parameter blocks for male outputs. The complete auditable both-sex rerun is `04_rerun_IN_comparison.R`; use it for the manuscript PDFs and numerical comparison.
