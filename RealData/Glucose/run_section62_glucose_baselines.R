#!/usr/bin/env Rscript

# RStudio-friendly Section 6.2 baseline runner.
# Open this file in RStudio and click Source to fit the LQD and Tangent
# baselines for the glucose example using the same seed, split, and
# preprocessing objects as the seed-20 proposed-method analysis.

get_this_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE), mustWork = TRUE))
  }
  source_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(source_path)) return(normalizePath(source_path, mustWork = TRUE))
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    path <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(path)) return(normalizePath(path, mustWork = TRUE))
  }
  stop("Run this file with Rscript or Source.")
}

script_path <- get_this_script_path()
script_dir <- dirname(script_path)

source(file.path(script_dir, "transport_fpca_modular.R"))
source(file.path(script_dir, "section62_glucose_helpers.R"))
source(file.path(script_dir, "baseline_fpca_methods.R"))

seed <- 20L
support_grid <- seq(50, 160, length.out = 201)
prob_grid <- seq(0, 1, length.out = 201)

Ku <- 10L
basis_type <- "bspline"
bs_degree <- 3L
methodXi <- "CE"
smooth_deriv <- FALSE
spar <- NULL

proposed_results_path <- file.path(
  script_dir,
  "output",
  sprintf("section62_glucose_seed_%s", seed),
  sprintf("section62_results_seed_%s.rds", seed)
)

if (!file.exists(proposed_results_path)) {
  stop(
    paste(
      "Saved proposed-method results not found at",
      proposed_results_path,
      "\nRun run_section62_glucose.R first."
    )
  )
}

proposed_results <- readRDS(proposed_results_path)
if (!identical(as.integer(proposed_results$seed), seed)) {
  stop(sprintf("Expected proposed results with seed %s, found seed %s.", seed, proposed_results$seed))
}

time_grid <- as.numeric(proposed_results$fit$workGrid)
data_path <- Sys.getenv("GLUCOSE_DATA_FILE", unset = file.path(script_dir, "Data", "data_1.csv"))
output_dir <- file.path(script_dir, "output", sprintf("section62_glucose_seed_%s", seed))

loaded <- glucose_load_data(data_path = data_path, day_window = proposed_results$config$day_window)
day_samples <- glucose_make_day_samples(loaded$data, loaded$subject_ids, loaded$keep_days)
subject_quantiles <- glucose_subject_quantiles(
  day_samples = day_samples,
  support_grid = support_grid,
  prob_grid = prob_grid,
  spline_spar = proposed_results$config$spline_spar
)

split_check <- glucose_make_split(
  keep_days = loaded$keep_days,
  seed = seed,
  holdout_days = proposed_results$config$holdout_days,
  min_days_for_holdout = proposed_results$config$min_days_for_holdout
)

if (!identical(split_check$train_days, proposed_results$split$train_days) ||
    !identical(split_check$test_days, proposed_results$split$test_days)) {
  stop("Recomputed split does not match the saved proposed-method split.")
}

transports <- glucose_make_transports(
  subject_quantiles = subject_quantiles,
  keep_days = loaded$keep_days,
  mean_cdf = proposed_results$mean_fit$mean_cdf,
  time_grid = time_grid,
  support_grid = support_grid
)

train_quantiles <- glucose_pick_rows(subject_quantiles, proposed_results$split$train_days)
train_transports <- glucose_pick_rows(transports$orig, proposed_results$split$train_days)
observed_holdouts <- glucose_pick_rows(transports$orig, proposed_results$split$test_days)

fpca_optns <- list(
  nRegGrid = length(time_grid),
  outPercent = c(0, 1),
  error = FALSE
)

lqd_fit <- fit_LQD_fpca(
  Q_list = train_quantiles,
  t_list = proposed_results$split$train_days,
  ugrid = prob_grid,
  Ku = Ku,
  basis_type = basis_type,
  bs_degree = bs_degree,
  t_grid = time_grid,
  methodXi = methodXi,
  fpca_optns = fpca_optns,
  smooth_deriv = smooth_deriv,
  spar = spar,
  lower = min(support_grid),
  upper = max(support_grid),
  project = TRUE,
  verbose = FALSE
)

lqd_transport_grid <- glucose_quantile_predictions_to_transports(
  Q_hat_list = lqd_fit$Q_hat,
  mean_cdf = proposed_results$mean_fit$mean_cdf,
  time_grid = time_grid,
  support_grid = support_grid,
  prob_grid = prob_grid
)

tangent_fit <- fit_tangent_wass_fpca(
  T_list = train_transports,
  t_list = proposed_results$split$train_days,
  xgrid = support_grid,
  Ku = Ku,
  basis_type = basis_type,
  bs_degree = bs_degree,
  t_grid = time_grid,
  methodXi = methodXi,
  fpca_optns = fpca_optns,
  lower = min(support_grid),
  upper = max(support_grid),
  project = TRUE,
  verbose = FALSE
)

lqd_holdouts <- glucose_select_prediction_times(
  curve_list = lqd_transport_grid,
  time_grid = time_grid,
  target_times = proposed_results$split$test_days
)
tangent_holdouts <- glucose_select_prediction_times(
  curve_list = tangent_fit$T_hat,
  time_grid = time_grid,
  target_times = proposed_results$split$test_days
)

lqd_errors <- glucose_holdout_error_table(
  pred_orig = lqd_holdouts,
  obs_orig = observed_holdouts,
  subject_ids = loaded$subject_ids,
  test_days = proposed_results$split$test_days,
  support_grid = support_grid
)
tangent_errors <- glucose_holdout_error_table(
  pred_orig = tangent_holdouts,
  obs_orig = observed_holdouts,
  subject_ids = loaded$subject_ids,
  test_days = proposed_results$split$test_days,
  support_grid = support_grid
)

baseline_results_path <- file.path(output_dir, sprintf("baseline_results_seed_%s.rds", seed))
lqd_errors_path <- file.path(output_dir, sprintf("baseline_lqd_errors_seed_%s.csv", seed))
tangent_errors_path <- file.path(output_dir, sprintf("baseline_tangent_errors_seed_%s.csv", seed))

saveRDS(
  list(
    seed = seed,
    config = list(
      Ku = Ku,
      basis_type = basis_type,
      bs_degree = bs_degree,
      methodXi = methodXi,
      smooth_deriv = smooth_deriv,
      spar = spar
    ),
    support_grid = support_grid,
    prob_grid = prob_grid,
    time_grid = time_grid,
    split = proposed_results$split,
    selected_summary = proposed_results$selected_summary,
    observed_holdouts = observed_holdouts,
    lqd_fit = lqd_fit,
    tangent_fit = tangent_fit,
    lqd_transport_grid = lqd_transport_grid,
    tangent_transport_grid = tangent_fit$T_hat,
    lqd_holdouts = lqd_holdouts,
    tangent_holdouts = tangent_holdouts,
    lqd_errors = lqd_errors,
    tangent_errors = tangent_errors,
    files = list(
      lqd_errors = lqd_errors_path,
      tangent_errors = tangent_errors_path
    )
  ),
  baseline_results_path
)

utils::write.csv(lqd_errors, lqd_errors_path, row.names = FALSE)
utils::write.csv(tangent_errors, tangent_errors_path, row.names = FALSE)

cat(sprintf("Finished Section 6.2 baseline analysis for seed %s\n", seed))
cat(sprintf("Saved baseline results: %s\n", baseline_results_path))
cat(sprintf("Saved LQD error table: %s\n", lqd_errors_path))
cat(sprintf("Saved Tangent error table: %s\n", tangent_errors_path))

selected_idx <- proposed_results$selected_summary$subject_index
selected_lqd <- lqd_errors[lqd_errors$subject_index %in% selected_idx, ]
selected_tangent <- tangent_errors[tangent_errors$subject_index %in% selected_idx, ]

cat("\nLQD errors for plotted subjects:\n")
print(selected_lqd[match(selected_idx, selected_lqd$subject_index), ], row.names = FALSE)

cat("\nTangent errors for plotted subjects:\n")
print(selected_tangent[match(selected_idx, selected_tangent$subject_index), ], row.names = FALSE)
