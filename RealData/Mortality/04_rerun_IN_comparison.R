# Run from any working directory:
# Rscript /path/to/Code/RealData/Mortality/04_rerun_IN_comparison.R [output_dir]
# Refit the two sexes with direct-U IN scores and a matched direct-U CE control.
# Original data and archived figures are never overwritten.

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) stop("Run this driver with Rscript.")
root <- dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg), fixed = TRUE), mustWork = TRUE))
source(file.path(root, "mortality_repro_functions.R"))
input_data_files <- mortality_input_files(root)
assert_packages(c("fdapace", "ggplot2"))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else file.path(root, "output", "estimator13_IN")
out <- ensure_dir(out)
if (file.exists(file.path(out, "run_complete.txt"))) {
  stop("Completed output already exists; supply a new output directory to rerun.")
}

config <- list(
  latent = FALSE, methods = c("IN", "CE"), kappa = 500,
  project_baselines = FALSE, u_cap = 0.999999, z_cap = NULL,
  include_mean = FALSE, components = 3L, score_method = "xi_est",
  legacy_baseline_postprocess = TRUE, project_reconstructed_rows = FALSE,
  fpca_options = list(methodMuCovEst = "smooth", usergrid = TRUE, methodSelectK = 3),
  country_codes_displayed = c("BLR", "ESP"),
  old_descriptive = list(latent = TRUE, kappa = 2000, method_xi = "IN"),
  metric_domain = "age 20--100; calendar years 1983--2018",
  integration = "Trapezoidal integration in age and calendar time; normalized metrics divide by domain lengths."
)
saveRDS(config, file.path(out, "configuration.rds"))
writeLines(capture.output(dput(config)), file.path(out, "configuration.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
inputs <- c(input_data_files,
            file.path(root, c("04_rerun_IN_comparison.R", "mortality_repro_functions.R",
                             "transport_fpca_modular.R")))
write.csv(data.frame(file = inputs, md5 = unname(tools::md5sum(inputs))),
          file.path(out, "input_code_md5.csv"), row.names = FALSE)

integrate_surface <- function(z, xgrid, times) {
  trapz_vec(times, apply(z, 1, function(v) trapz_vec(xgrid, v)))
}

reconstruct_for_comparison <- function(fit, context) {
  # Preserve the baseline postprocessing and mean omission of the archived
  # trajectory script; the paired fit differs only in methodXi.
  upper <- max(context$xgrid)
  fit$T0_plus <- t(apply(fit$T0_plus, 1, legacy_monotone_upper_clamp, upper_bound = upper))
  fit$T0_minus <- t(apply(fit$T0_minus, 1, legacy_monotone_upper_clamp, upper_bound = upper))
  scalar <- reconstruct_scalar_trajectories(fit, include_mean = FALSE, score_method = "xi_est")
  list(fit = fit, scalar = scalar,
       transports = reconstruct_transports_from_scalar(fit, scalar, project_rows = FALSE))
}

comparison_metrics <- function(context, IN, CE) {
  age_range <- diff(range(context$xgrid))
  time_range <- diff(range(context$years))
  domain_size <- age_range * time_range
  rows <- lapply(seq_along(context$country_codes), function(i) {
    observed <- context$transport_list[[i]]
    deviation <- sweep(observed, 2, context$xgrid, "-")
    delta <- IN$transports[[i]] - CE$transports[[i]]
    integrated <- function(z) integrate_surface(z, context$xgrid, context$years)
    iae <- integrated(abs(delta))
    ise <- integrated(delta^2)
    observed_abs <- integrated(abs(deviation))
    observed_sq <- integrated(deviation^2)
    in_sq <- integrated((IN$transports[[i]] - observed)^2)
    ce_sq <- integrated((CE$transports[[i]] - observed)^2)
    data.frame(
      sex = context$sex, country_code = context$country_codes[i],
      country = context$country_labels[i],
      integrated_abs_IN_minus_CE = iae, integrated_sq_IN_minus_CE = ise,
      mean_abs_IN_minus_CE_years = iae / domain_size,
      rms_IN_minus_CE_years = sqrt(ise / domain_size),
      max_abs_IN_minus_CE_years = max(abs(delta)),
      observed_mean_abs_deviation_years = observed_abs / domain_size,
      observed_rms_deviation_years = sqrt(observed_sq / domain_size),
      abs_difference_pct_observed_abs = 100 * iae / observed_abs,
      rms_difference_pct_observed_rms = 100 * sqrt(ise / observed_sq),
      IN_in_sample_RMSE_years = sqrt(in_sq / domain_size),
      CE_in_sample_RMSE_years = sqrt(ce_sq / domain_size),
      IN_in_sample_RMSE_pct_observed_rms = 100 * sqrt(in_sq / observed_sq),
      CE_in_sample_RMSE_pct_observed_rms = 100 * sqrt(ce_sq / observed_sq),
      IN_minus_CE_in_sample_MSE_pct_CE = 100 * (in_sq / ce_sq - 1),
      scalar_max_abs_difference = max(abs(IN$scalar[i, ] - CE$scalar[i, ])),
      scalar_sign_agreement = mean(sign(IN$scalar[i, ]) == sign(CE$scalar[i, ]))
    )
  })
  do.call(rbind, rows)
}

save_trajectory_pdf <- function(context, transports, country_code, filename) {
  i <- match(country_code, context$country_codes)
  stopifnot(!is.na(i))
  # Only the surface mesh is thinned for a compact vector PDF; metrics use all
  # 401 age-grid points and all 36 annual observations.
  ages_idx <- seq(1, length(context$xgrid), by = 4)
  ages <- context$xgrid[ages_idx]
  fitted <- sweep(transports[[i]], 2, context$xgrid, "-")[, ages_idx]
  observed <- sweep(context$transport_list[[i]], 2, context$xgrid, "-")[, ages_idx]
  grDevices::pdf(filename, width = 6.4, height = 5.7, family = "Helvetica", useDingbats = FALSE)
  on.exit(grDevices::dev.off())
  graphics::layout(matrix(1:2, ncol = 1), heights = c(0.92, 0.08))
  graphics::par(mar = c(2, 1.3, 2.4, 1.8), mgp = c(1.8, 0.4, 0), cex = 0.95)
  projection <- graphics::persp(
    x = context$years, y = ages, z = fitted,
    zlim = range(c(fitted, observed, 0)), theta = -40, phi = 23,
    expand = 0.7, ticktype = "detailed", nticks = 4,
    col = "#8BAAF0", border = NA, shade = 0.2,
    xlab = "Calendar year", ylab = "Age", zlab = "Transport - identity",
    main = paste0(country_code_to_label(country_code), " (", country_code, ")")
  )
  for (j in seq_along(context$years)) {
    projected <- grDevices::trans3d(rep(context$years[j], length(ages)), ages,
                                  observed[j, ], pmat = projection)
    graphics::lines(projected, col = grDevices::adjustcolor("#CB2525", alpha.f = 0.75), lwd = 0.75)
  }
  graphics::par(mar = rep(0, 4))
  graphics::plot.new()
  graphics::legend("center", legend = c("IN fitted surface", "Observed transports"),
                   fill = c("#8BAAF0", NA), border = c(NA, NA),
                   lty = c(NA, 1), col = c(NA, "#CB2525"), bty = "n", horiz = TRUE, cex = 1.5)
}

all_metrics <- list()
eigenvalues <- list()
full_spectra <- list()
diagnostics <- list()
old_comparison <- list()
invariance <- list()
signs <- list()
female_order <- c("JPN", "SWE", "ITA", "ESP", "NOR", "ISL", "CHE", "NLD", "GRC",
                  "AUS", "FRATNP", "AUT", "LUX", "CAN", "FIN", "SVN", "GBRTENW",
                  "PRT", "IRL", "BEL", "NZL_NP", "DNK", "CZE", "TWN", "POL", "SVK",
                  "USA", "EST", "BGR", "HUN", "LTU", "BLR", "LVA")

for (sex in c("Female", "Male")) {
  cat("Running", sex, "\n")
  context <- create_mortality_repro_context(root, sex)
  stopifnot(length(context$country_codes) == 33L, identical(context$years, 1983:2018),
            length(context$xgrid) == 401L, all(is.finite(unlist(context$transport_list))))
  fit <- lapply(config$methods, function(method) {
    fit_mortality_transport_fpca(context, kappa = config$kappa, latent = FALSE,
                                method_xi = method, project_baselines = FALSE,
                                u_cap = config$u_cap, z_cap = NULL,
                                fpca_options = config$fpca_options)
  })
  names(fit) <- config$methods
  stopifnot(all(vapply(fit, function(f) f$fpca_fit$selectK == 3L, logical(1))))
  check <- c(mu = max(abs(fit$IN$fpca_fit$mu - fit$CE$fpca_fit$mu)),
             covariance = max(abs(fit$IN$fpca_fit$smoothedCov - fit$CE$fpca_fit$smoothedCov)),
             eigenvalues = max(abs(fit$IN$fpca_fit$lambda - fit$CE$fpca_fit$lambda)),
             eigenfunctions = max(abs(fit$IN$fpca_fit$phi - fit$CE$fpca_fit$phi)))
  stopifnot(all(check < 1e-12))
  invariance[[sex]] <- data.frame(sex = sex, quantity = names(check), max_abs_difference = unname(check))
  recon <- lapply(fit, reconstruct_for_comparison, context = context)
  all_metrics[[sex]] <- comparison_metrics(context, recon$IN, recon$CE)
  f <- fit$IN$fpca_fit
  # Recompute the positive covariance spectrum exactly as fdapace does (up to
  # maxK), so the denominator behind the returned cumFVE can be audited.
  spectrum <- eigen(f$smoothedCov, symmetric = TRUE, only.values = TRUE)$values
  spectrum <- head(spectrum[spectrum >= 0], f$optns$maxK) * diff(f$workGrid)[1]
  stopifnot(max(abs(cumsum(spectrum)[1:3] / sum(spectrum) - f$cumFVE[1:3])) < 1e-10)
  eigenvalues[[sex]] <- data.frame(sex = sex, component = 1:3, eigenvalue = f$lambda,
                                  individual_FVE_pct = 100 * f$lambda / sum(spectrum),
                                  cumulative_FVE_pct = 100 * f$cumFVE[1:3],
                                  total_positive_eigenvalue_sum = sum(spectrum))
  full_spectra[[sex]] <- data.frame(sex = sex, component = seq_along(spectrum), eigenvalue = spectrum)
  d <- fit$IN$diagnostics
  diagnostics[[sex]] <- data.frame(sex = sex, subjects = 33, observations = 1188,
                                   clipped_observations = sum(abs(unlist(fit$IN$U_hat_raw)) > config$u_cap),
                                   clipped_fraction = mean(abs(unlist(fit$IN$U_hat_raw)) > config$u_cap),
                                   max_abs_U_raw = d$max_abs_U_raw,
                                   estimated_mean_min = min(f$mu), estimated_mean_max = max(f$mu),
                                   data_type = f$optns$dataType, kernel = f$optns$kernel,
                                   bandwidth_mean = f$bwMu, bandwidth_covariance = f$bwCov)
  old <- fit_mortality_transport_fpca(context, kappa = 2000, latent = TRUE, method_xi = "IN")
  phi_corr <- vapply(1:3, function(k) {
    a <- f$phi[, k]; b <- old$fpca_fit$phi[, k]
    abs(trapz_vec(f$workGrid, a * b) /
          sqrt(trapz_vec(f$workGrid, a^2) * trapz_vec(f$workGrid, b^2)))
  }, numeric(1))
  old_comparison[[sex]] <- data.frame(sex = sex, component = 1:3,
                                     absolute_L2_inner_product = phi_corr,
                                     new_direct_U_cumulative_FVE_pct = 100 * f$cumFVE[1:3],
                                     old_latent_cumulative_FVE_pct = 100 * old$fpca_fit$cumFVE[1:3],
                                     old_latent_mean_min = min(old$fpca_fit$mu),
                                     old_latent_mean_max = max(old$fpca_fit$mu))
  sex_tag <- sex_suffix(sex)
  mean_plot <- build_mean_eigen_plot(fit$IN, use_legacy_fve_labels = FALSE,
                                    eigen_sign_flip = c(-1, 1, -1),
                                    show_axis_lines = FALSE,
                                    legend_title_size = 11, legend_text_size = 11)
  save_ggplot_pdf(mean_plot, file.path(out, paste0("mean_IN_", sex_tag, ".pdf")),
                  width = 7.0, height = 5.5, use_cairo = FALSE)
  heatmap <- run_country_sign_figure(context, output_dir = out,
                                    output_filename = paste0("countries_", sex_tag, ".pdf"),
                                    ordered_country_codes = if (sex == "Female") female_order else NULL,
                                    use_cairo_pdf = FALSE, display_plot = FALSE)
  signs[[sex]] <- data.frame(sex = sex,
                             country_code = rep(context$country_codes, each = length(context$years)),
                             year = rep(context$years, times = length(context$country_codes)),
                             signed_mass = as.vector(t(context$signed_mass_matrix)))
  if (sex == "Female") {
    for (country in config$country_codes_displayed) {
      save_trajectory_pdf(context, recon$IN$transports, country,
                          file.path(out, paste0(tolower(country), "_IN_F.pdf")))
    }
  }
  saveRDS(list(configuration = config, countries = context$country_codes,
               years = context$years, xgrid = context$xgrid,
               observed_transports = context$transport_list, signed_mass = context$signed_mass_matrix,
               IN = recon$IN, CE = recon$CE, old_latent_descriptive_fit = old),
          file.path(out, paste0("mortality_", sex_tag, "_fits.rds")))
}

metrics <- do.call(rbind, all_metrics)
write.csv(metrics, file.path(out, "comparison_by_country.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(metrics, metrics$sex), function(d) {
  data.frame(sex = d$sex[1], countries = nrow(d),
             mean_abs_difference_years = mean(d$mean_abs_IN_minus_CE_years),
             pooled_rms_difference_years = sqrt(mean(d$rms_IN_minus_CE_years^2)),
             maximum_pointwise_difference_years = max(d$max_abs_IN_minus_CE_years),
             pooled_abs_difference_pct_observed = 100 * sum(d$mean_abs_IN_minus_CE_years) /
               sum(d$observed_mean_abs_deviation_years),
             pooled_rms_difference_pct_observed = 100 * sqrt(sum(d$rms_IN_minus_CE_years^2) /
               sum(d$observed_rms_deviation_years^2)),
             max_country_rms_difference_pct_observed = max(d$rms_difference_pct_observed_rms),
             pooled_IN_in_sample_RMSE_years = sqrt(mean(d$IN_in_sample_RMSE_years^2)),
             pooled_CE_in_sample_RMSE_years = sqrt(mean(d$CE_in_sample_RMSE_years^2)),
             countries_with_lower_IN_in_sample_MSE = sum(d$IN_in_sample_RMSE_years < d$CE_in_sample_RMSE_years),
             scalar_sign_agreement = mean(d$scalar_sign_agreement))
}))
write.csv(summary, file.path(out, "comparison_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, eigenvalues), file.path(out, "fitted_eigenvalues_FVE.csv"), row.names = FALSE)
write.csv(do.call(rbind, full_spectra), file.path(out, "positive_covariance_spectrum.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostics), file.path(out, "fit_diagnostics.csv"), row.names = FALSE)
write.csv(do.call(rbind, old_comparison), file.path(out, "descriptive_old_latent_vs_direct_U.csv"), row.names = FALSE)
write.csv(do.call(rbind, invariance), file.path(out, "CE_IN_covariance_invariance.csv"), row.names = FALSE)
write.csv(do.call(rbind, signs), file.path(out, "empirical_signed_mass_by_country_year.csv"), row.names = FALSE)
write.csv(data.frame(code = c("BLR", "ESP"), country = c("Belarus", "Spain"),
                     original_html = c("blr_F.html", "bel_F.html"),
                     updated_pdf = c("blr_IN_F.pdf", "esp_IN_F.pdf")),
          file.path(out, "displayed_country_mapping.csv"), row.names = FALSE)
writeLines(c("All 66 country/sex matched comparisons completed.",
             "CE and IN means, covariance estimates, eigenvalues, and eigenfunctions agree exactly.",
             "FVE checked against the fitted positive covariance spectrum, before truncation to three PCs.",
             "Errors are in-sample transport discrepancies, not held-out prediction errors.",
             "The direct-U IN implementation uses trapezoidal integration of centered observations.",
             "Original baseline postprocessing, mean omission, and observation clipping are retained.",
             paste("Completed:", format(Sys.time(), tz = "UTC", usetz = TRUE))),
           file.path(out, "run_complete.txt"))
print(summary)
cat("Outputs:", out, "\n")
