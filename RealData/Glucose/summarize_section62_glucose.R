#!/usr/bin/env Rscript

# Audit and aggregate the saved Section 6.2 holdout predictions, without refitting.
# Run after run_section62_glucose.R and run_section62_glucose_baselines.R.
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
source(file.path(script_dir, "section62_glucose_helpers.R"))
seed <- 20L
output_dir <- file.path(script_dir, "output", sprintf("section62_glucose_seed_%s", seed))
proposed_path <- file.path(output_dir, sprintf("section62_results_seed_%s.rds", seed))
baseline_path <- file.path(output_dir, sprintf("baseline_results_seed_%s.rds", seed))
data_path <- Sys.getenv("GLUCOSE_DATA_FILE", unset = file.path(script_dir, "Data", "data_1.csv"))
if (!file.exists(proposed_path) || !file.exists(baseline_path)) {
  stop("Run run_section62_glucose.R and run_section62_glucose_baselines.R first.", call. = FALSE)
}
if (!file.exists(data_path)) {
  stop("Supply Data/data_1.csv or set GLUCOSE_DATA_FILE; see the application README.", call. = FALSE)
}
proposed <- readRDS(proposed_path)
baseline <- readRDS(baseline_path)

assert_equal <- function(x, y, label, tolerance = 1e-10) {
  if (!isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))) {
    stop(sprintf("Validation failed: %s", label))
  }
}

stopifnot(proposed$seed == seed, baseline$seed == seed,
          proposed$config$min_days_for_holdout == 6,
          proposed$config$holdout_days == 2,
          isTRUE(proposed$config$latent), isTRUE(proposed$fit$latent),
          identical(baseline$config$methodXi, "CE"))
proposed_score <- proposed$fit$fpca_fit$optns$methodXi
if (is.null(proposed_score)) proposed_score <- as.character(proposed$fit$call$methodXi)
stopifnot(identical(proposed_score, "CE"))
assert_equal(proposed$split, baseline$split, "common train/test split")
assert_equal(proposed$observed_holdouts, baseline$observed_holdouts, "common observed targets")
assert_equal(proposed$selected_summary, baseline$selected_summary, "common selected subjects")

# Recreate only the deterministic data preparation and split to audit eligibility.
loaded <- glucose_load_data(data_path, day_window = proposed$config$day_window)
assert_equal(loaded$subject_ids, proposed$subject_ids, "subject identifiers")
assert_equal(loaded$keep_days, proposed$keep_days, "eligible study days")
split_check <- glucose_make_split(loaded$keep_days, seed, holdout_days = 2,
                                 min_days_for_holdout = 6)
assert_equal(split_check, proposed$split, "seed-20 split reproduced from supplied data")
eligible <- unname(which(lengths(proposed$keep_days) >= 6))
stopifnot(length(eligible) > 1,
          identical(unname(which(lengths(proposed$split$test_days) > 0)), eligible),
          all(lengths(proposed$split$test_days)[eligible] == 2))
for (i in seq_along(proposed$subject_ids)) {
  train <- proposed$split$train_days[[i]]
  test <- proposed$split$test_days[[i]]
  stopifnot(!anyDuplicated(c(train, test)))
  assert_equal(sort(c(train, test)), sort(proposed$keep_days[[i]]),
               sprintf("train/test partition for subject %s", i))
}

# Saved curves are on the original support, whereas proposed$fit$xgrid is scaled.
x <- as.numeric(baseline$support_grid)
stopifnot(length(x) > 1, all(is.finite(x)), all(diff(x) > 0),
          min(x) == 50, max(x) == 160)
assert_equal(proposed$fit$xgrid, (x - min(x)) / diff(range(x)), "scaled FPCA grid")
predictions <- list(Proposed = proposed$predicted_holdouts,
                    LQD = baseline$lqd_holdouts, Tangent = baseline$tangent_holdouts)
observed <- proposed$observed_holdouts
selected <- as.integer(proposed$selected_summary$subject_index)
stopifnot(length(selected) == 2, !anyDuplicated(selected), all(selected %in% eligible))
rows <- list()
for (method in names(predictions)) {
  stopifnot(length(predictions[[method]]) == length(proposed$subject_ids))
  for (i in eligible) {
    days <- proposed$split$test_days[[i]]
    obs <- observed[[i]]
    pred <- predictions[[method]][[i]]
    stopifnot(is.matrix(obs), is.matrix(pred),
              identical(dim(obs), c(2L, length(x))), identical(dim(pred), dim(obs)),
              identical(rownames(obs), as.character(days)),
              identical(rownames(pred), as.character(days)),
              all(is.finite(obs)), all(is.finite(pred)))
    for (j in seq_along(days)) {
      delta <- abs(pred[j, ] - obs[j, ])
      error <- sum(diff(x) * (head(delta, -1) + tail(delta, -1)) / 2)
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, subject_index = i, subject_id = proposed$subject_ids[i],
        holdout_day = days[j], selected_for_figure = i %in% selected,
        integrated_absolute_error = error, squared_integrated_absolute_error = error^2,
        stringsAsFactors = FALSE
      )
    }
  }
}
per_day <- do.call(rbind, rows)
stopifnot(nrow(per_day) == 3L * 2L * length(eligible),
          all(is.finite(per_day$integrated_absolute_error)),
          all(is.finite(per_day$squared_integrated_absolute_error)))

per_subject <- do.call(rbind, lapply(names(predictions), function(method) {
  do.call(rbind, lapply(eligible, function(i) {
    z <- per_day[per_day$method == method & per_day$subject_index == i, ]
    data.frame(method = method, subject_index = i, subject_id = z$subject_id[1],
               n_holdout_days = nrow(z), selected_for_figure = i %in% selected,
               mean_integrated_absolute_error = mean(z$integrated_absolute_error),
               mean_squared_integrated_absolute_error = mean(z$squared_integrated_absolute_error),
               stringsAsFactors = FALSE)
  }))
}))

# Check full-precision means and ranking; rounded comma-separated errors are unused.
saved_errors <- list(Proposed = proposed$error_table,
                     LQD = baseline$lqd_errors, Tangent = baseline$tangent_errors)
for (method in names(predictions)) {
  z <- per_subject[per_subject$method == method, ]
  saved <- saved_errors[[method]]
  stopifnot(nrow(saved) == length(eligible), !anyDuplicated(saved$subject_index))
  assert_equal(sort(as.integer(saved$subject_index)), eligible, paste(method, "complete subject set"))
  saved_order <- match(z$subject_index, saved$subject_index)
  assert_equal(saved$subject_id[saved_order], z$subject_id, paste(method, "saved subject IDs"))
  assert_equal(saved$mean_holdout_w1[saved_order], z$mean_integrated_absolute_error,
               paste(method, "saved full-precision errors"))
  rank_order <- order(z$mean_integrated_absolute_error)
  assert_equal(saved$subject_index, z$subject_index[rank_order], paste(method, "saved ranking"))
  if (method == "Proposed") {
    assert_equal(selected, head(z$subject_index[rank_order], 2), "two minimum-error subjects selected")
  }
}

aggregate_rows <- list()
for (scope in c("all_eligible", "selected_for_figure")) {
  for (method in names(predictions)) {
    z <- per_subject[per_subject$method == method, ]
    if (scope == "selected_for_figure") z <- z[z$selected_for_figure, ]
    a <- z$mean_integrated_absolute_error
    b <- z$mean_squared_integrated_absolute_error
    aggregate_rows[[length(aggregate_rows) + 1L]] <- data.frame(
      scope = scope, method = method, n_subjects = nrow(z), n_holdout_days = sum(z$n_holdout_days),
      mean_integrated_absolute_error = mean(a), subject_se_integrated_absolute_error = sd(a) / sqrt(length(a)),
      median_subject_integrated_absolute_error = median(a),
      mean_squared_integrated_absolute_error = mean(b),
      subject_se_squared_integrated_absolute_error = sd(b) / sqrt(length(b)),
      median_subject_squared_integrated_absolute_error = median(b), stringsAsFactors = FALSE
    )
  }
}
aggregate_results <- do.call(rbind, aggregate_rows)
utils::write.csv(per_day, file.path(output_dir, "holdout_per_day.csv"), row.names = FALSE)
utils::write.csv(per_subject, file.path(output_dir, "holdout_per_subject.csv"), row.names = FALSE)
utils::write.csv(aggregate_results, file.path(output_dir, "holdout_aggregate.csv"), row.names = FALSE)

source_files <- file.path(script_dir, c("run_section62_glucose.R", "run_section62_glucose_baselines.R",
                                      "section62_glucose_helpers.R", "transport_fpca_modular.R",
                                      "baseline_fpca_methods.R", "summarize_section62_glucose.R"))
audit <- c(
  "Section 6.2 glucose holdout audit", sprintf("Created: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("Session timezone: %s; TZ environment: %s", Sys.timezone(), Sys.getenv("TZ")),
  sprintf("Timestamp parsing timezone: %s; the data and split were reproduced.", Sys.getenv("GLUCOSE_TIMEZONE", "America/Chicago")),
  sprintf("Seed: %s; total subjects: %s; eligible subjects: %s; held-out days: %s",
          seed, length(proposed$subject_ids), length(eligible), 2 * length(eligible)),
  "Eligibility: at least six observed study days; exactly two days sampled without replacement per eligible subject.",
  "Proposed procedure: latent=TRUE, methodXi=CE (estimator (14)). LQD and Tangent use CE scores.",
  "Validation: identical splits and observed targets; complete finite predictions for every eligible subject/day;",
  "saved full-precision errors and rankings verified; selected subjects are the two lowest proposed mean errors.",
  "Error on day t: e_it = integral from 50 to 160 of |T_hat_i(t,x) - T_observed_i(t,x)| dx (trapezoidal rule).",
  "This is the raw integrated absolute transport error on the original glucose support; no normalization is applied.",
  "Squared error is e_it^2 (square AFTER integration), not the integral of squared pointwise differences.",
  "Subject errors average the two held-out days. Reported method means weight all eligible subjects equally.",
  "Subject-level SE = sample SD of subject means / sqrt(number of subjects); it is conditional on this split.",
  "The selected-subject summary is illustrative and is not an unbiased assessment of predictive performance.",
  "", "Selected subjects:", capture.output(print(proposed$selected_summary, row.names = FALSE)),
  "", "Aggregate results:", capture.output(print(aggregate_results, row.names = FALSE, digits = 12)),
  "", "Proposed configuration:", capture.output(dput(proposed$config)),
  "", "Baseline configuration:", capture.output(dput(baseline$config)),
  "", "Input/result/source MD5 hashes:",
  capture.output(print(tools::md5sum(c(data_path, proposed_path, baseline_path, source_files)))),
  "", "Session information:", capture.output(sessionInfo())
)
writeLines(audit, file.path(output_dir, "holdout_aggregation_audit.txt"))
print(aggregate_results, row.names = FALSE, digits = 8)
cat(sprintf("\nValidated all %s eligible subjects (%s held-out days per method).\n", length(eligible), 2 * length(eligible)))
cat(sprintf("Saved aggregate tables and audit in %s\n", output_dir))
