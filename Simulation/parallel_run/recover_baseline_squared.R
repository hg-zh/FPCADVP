#!/usr/bin/env Rscript

# Recover squared-W1 IMSE from archived pointwise mean squared-error curves.
# Usage: Rscript recover_baseline_squared.R input_dir=PATH output_dir=NEW_PATH
# Sourcing this file defines the recovery function without writing any files.
recover_baseline_squared <- function(input_dir, output_dir) {
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  output_dir <- path.expand(output_dir)
  output_names <- c("raw_results.rds", "recovery_validation.csv",
                    "source_checksums.csv", "recovery_note.txt")
  if (any(file.exists(file.path(output_dir, output_names)))) {
    stop("Recovery outputs already exist; choose a new output_dir.")
  }
  raw_path <- file.path(input_dir, "raw_results.rds")
  raw <- readRDS(raw_path)
  keys <- c("rep", "n", "N", "m")
  required <- c(keys, "ok", "imse_lqd", "imse_tangent")
  if (!is.data.frame(raw) || !all(required %in% names(raw))) {
    stop("Raw baseline results lack the required columns.")
  }
  if (any(c("imse_sq_lqd", "imse_sq_tangent") %in% names(raw))) {
    stop("Squared-error columns already exist; recovery must use the original archive.")
  }
  if (nrow(raw) != 4000L || !is.logical(raw$ok) || anyNA(raw$ok) || !all(raw$ok)) {
    stop("Expected exactly 4,000 successful baseline records.")
  }
  numeric_cols <- c(keys, "imse_lqd", "imse_tangent")
  if (!all(vapply(raw[numeric_cols], function(x) {
    is.numeric(x) && all(is.finite(x)) && all(x >= 0)
  }, logical(1)))) stop("Keys and legacy errors must be finite, nonnegative numbers.")
  key <- function(x) do.call(paste, c(x[keys], sep = "/"))
  expected <- expand.grid(rep = seq_len(200L), n = c(20, 50, 100, 200, 400),
                          N = c(3, 5, 10, 20), m = 200)
  if (anyDuplicated(key(raw)) || !setequal(key(raw), key(expected))) {
    stop("Raw keys must cover every planned (rep,n,N,m) combination exactly once.")
  }
  expected_names <- sprintf("MSE_curve_rep%03d_n%d_N%d_m%d.csv",
                            raw$rep, raw$n, raw$N, raw$m)
  curve_files <- list()
  validation <- list()
  for (method in c("LQD", "Tangent")) {
    curve_dir <- file.path(input_dir, "mse_curves", method)
    found <- list.files(curve_dir, pattern = "[.]csv$")
    if (length(found) != 4000L || !setequal(found, expected_names)) {
      stop(method, ": curve files do not cover all 4,000 planned records exactly.")
    }
    files <- file.path(curve_dir, expected_names)
    recovered <- numeric(nrow(raw))
    checks <- vector("list", nrow(raw))
    for (i in seq_len(nrow(raw))) {
      z <- utils::read.csv(files[i], stringsAsFactors = FALSE)
      if (!all(c("t", "mse") %in% names(z)) || nrow(z) != 51L ||
          !is.numeric(z$t) || !is.numeric(z$mse) ||
          any(!is.finite(z$t)) || any(!is.finite(z$mse)) || any(z$mse < 0) ||
          any(abs(z$t - seq(0, 1, length.out = 51L)) > 1e-12)) {
        stop("Invalid pointwise squared-error curve: ", files[i])
      }
      # The saved mse is mean_i W1_i(t)^2. Linear trapezoidal integration
      # commutes with averaging over subjects; do not square an averaged IW1.
      recovered[i] <- sum(diff(z$t) * (head(z$mse, -1L) + tail(z$mse, -1L)) / 2)
      if (!is.finite(recovered[i])) stop("Nonfinite recovered IMSE: ", files[i])
      checks[[i]] <- data.frame(raw[i, keys, drop = FALSE], method = method,
                               time_points = nrow(z), time_start = z$t[1L],
                               time_end = tail(z$t, 1L),
                               min_pointwise_mse = min(z$mse),
                               max_pointwise_mse = max(z$mse),
                               imse_sq = recovered[i], valid = TRUE)
    }
    raw[[paste0("imse_sq_", tolower(method))]] <- recovered
    curve_files[[method]] <- files
    validation[[method]] <- do.call(rbind, checks)
  }
  source_files <- c(raw_path, unlist(curve_files, use.names = FALSE))
  checksums <- data.frame(
    source = substring(source_files, nchar(input_dir) + 2L),
    md5 = unname(tools::md5sum(source_files)), stringsAsFactors = FALSE)
  if (anyNA(checksums$md5)) stop("Could not checksum every source file.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output directory: ", output_dir)
  # Recheck immediately before writing to protect concurrently created outputs.
  if (any(file.exists(file.path(output_dir, output_names)))) {
    stop("Recovery outputs appeared during validation; no files were overwritten.")
  }
  saveRDS(raw, file.path(output_dir, "raw_results.rds"))
  utils::write.csv(do.call(rbind, validation),
                   file.path(output_dir, "recovery_validation.csv"), row.names = FALSE)
  utils::write.csv(checksums, file.path(output_dir, "source_checksums.csv"), row.names = FALSE)
  writeLines(c(
    "Squared-W1 baseline IMSE recovered from archived pointwise MSE curves.",
    paste("Source directory:", input_dir),
    "All 4,000 successful records and all 8,000 expected curves passed validation.",
    "Grid: 200 repetitions; n=20,50,100,200,400; N=3,5,10,20; m=200.",
    "Each curve has 51 finite, nonnegative mean_i W1_i(t)^2 values on t=0:0.02:1.",
    "New columns imse_sq_lqd and imse_sq_tangent integrate these values by trapezoidal quadrature.",
    "All original raw columns, including unsquared imse_lqd/imse_tangent, are preserved.",
    "No model fits were rerun; no unsquared aggregate was squared.",
    "Source checksums identify the exact archives. This recovery verifies the metric and completeness;",
    "the archived data-generating settings must be verified separately from their run provenance."
  ), file.path(output_dir, "recovery_note.txt"))
  message("Recovered squared IMSE for 4,000 records: ", normalizePath(output_dir))
  invisible(raw)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L || any(!grepl("=", args, fixed = TRUE))) {
    stop("Usage: Rscript recover_baseline_squared.R input_dir=PATH output_dir=NEW_PATH")
  }
  names_arg <- sub("=.*$", "", args)
  if (anyDuplicated(names_arg) || !setequal(names_arg, c("input_dir", "output_dir"))) {
    stop("Supply exactly input_dir=PATH and output_dir=NEW_PATH.")
  }
  values <- setNames(sub("^[^=]*=", "", args), names_arg)
  if (any(!nzchar(values))) stop("Input and output directories must be nonempty.")
  recover_baseline_squared(values[["input_dir"]], values[["output_dir"]])
}
