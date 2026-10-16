#!/usr/bin/env Rscript

## Kernel/bandwidth sensitivity driver for the proposed transport-FPCA method.
##
## The scientific defaults implement Setting 1 in Section 5 of the manuscript.
## Every (rep, n, N, m) dataset is generated once.  The untouched implicit
## fdapace fit is run first and reused as the gauss x 1.00 reference.  All
## other fits use explicit kernel, userBwMu, and userBwCov options obtained by
## matching the standard kernels' second moments and then multiplying that
## dataset's returned Gaussian reference bandwidths by the sensitivity factor.

options(stringsAsFactors = FALSE)

suppressMessages(library(parallel))

`%||%` <- function(a, b) if (!is.null(a)) a else b

get_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) > 0L) {
    script_path <- gsub("~+~", " ", sub("^--file=", "", file_arg[1L]), fixed = TRUE)
    return(dirname(normalizePath(script_path, mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

parse_kv_args <- function(args) {
  ans <- list()
  if (length(args) == 0L) return(ans)
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) {
      stop("Every command-line argument must have key=value form; invalid argument: ", arg)
    }
    key <- sub("=.*$", "", arg)
    value <- sub("^[^=]*=", "", arg)
    if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", key)) {
      stop("Invalid command-line key: ", key)
    }
    if (!is.null(ans[[key]])) stop("Duplicate command-line key: ", key)
    ans[[key]] <- value
  }
  ans
}

parse_bool <- function(x, default = FALSE, name = "value") {
  if (is.null(x)) return(default)
  value <- tolower(trimws(as.character(x)[1L]))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(name, " must be one of TRUE/FALSE, T/F, 1/0, yes/no")
}

parse_int_scalar <- function(x, default, name, lower = -Inf) {
  if (is.null(x)) return(as.integer(default))
  value_num <- suppressWarnings(as.numeric(trimws(x)))
  value <- suppressWarnings(as.integer(value_num))
  if (!is.finite(value_num) || is.na(value) || value_num != value || value < lower) {
    stop(name, " must be an integer >= ", lower)
  }
  value
}

parse_num_scalar <- function(x, default, name, lower = -Inf, upper = Inf,
                             lower_open = FALSE, upper_open = FALSE) {
  if (is.null(x)) return(as.numeric(default))
  value <- suppressWarnings(as.numeric(trimws(x)))
  lower_bad <- if (lower_open) value <= lower else value < lower
  upper_bad <- if (upper_open) value >= upper else value > upper
  if (length(value) != 1L || !is.finite(value) || lower_bad || upper_bad) {
    stop(name, " is outside its allowed numeric range")
  }
  value
}

split_csv <- function(x) {
  value <- gsub("\\s+", "", paste(as.character(x), collapse = ""))
  if (!nzchar(value)) return(character(0))
  strsplit(value, ",", fixed = TRUE)[[1L]]
}

parse_int_vec <- function(x, default, name) {
  values <- if (is.null(x)) as.character(default) else split_csv(x)
  nums <- suppressWarnings(as.numeric(values))
  out <- suppressWarnings(as.integer(nums))
  if (length(out) < 1L || any(!is.finite(nums)) || any(is.na(out)) ||
      any(nums != out) || any(out <= 0L) || anyDuplicated(out)) {
    stop(name, " must be a comma-separated vector of distinct positive integers")
  }
  out
}

parse_num_vec <- function(x, default, name) {
  values <- if (is.null(x)) as.character(default) else split_csv(x)
  out <- suppressWarnings(as.numeric(values))
  if (length(out) < 1L || any(!is.finite(out)) || any(out <= 0) || anyDuplicated(out)) {
    stop(name, " must be a comma-separated vector of distinct positive numbers")
  }
  out
}

parse_char_vec <- function(x, default, name) {
  out <- if (is.null(x)) as.character(default) else split_csv(x)
  if (length(out) < 1L || any(!nzchar(out)) || anyDuplicated(out)) {
    stop(name, " must be a comma-separated vector of distinct nonempty values")
  }
  out
}

same_numeric_values <- function(x, y, tolerance = 1e-14) {
  length(x) == length(y) &&
    isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance,
                     check.attributes = FALSE))
}

close_numeric <- function(x, y, tolerance = 1e-10) {
  length(x) == 1L && length(y) == 1L && is.finite(x) && is.finite(y) &&
    abs(x - y) <= tolerance * max(1, abs(x), abs(y))
}

atomic_save_rds <- function(object, path) {
  dir_path <- dirname(path)
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) stop("Refusing to overwrite existing file: ", path)
  tmp <- tempfile(pattern = paste0(".", basename(path), ".tmp-"), tmpdir = dir_path)
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) stop("Atomic rename failed for: ", path)
  invisible(path)
}

atomic_write_csv <- function(object, path) {
  dir_path <- dirname(path)
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) stop("Refusing to overwrite existing file: ", path)
  tmp <- tempfile(pattern = paste0(".", basename(path), ".tmp-"), tmpdir = dir_path)
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  utils::write.csv(object, tmp, row.names = FALSE, na = "NA")
  if (!file.rename(tmp, path)) stop("Atomic rename failed for: ", path)
  invisible(path)
}

atomic_write_lines <- function(text, path) {
  dir_path <- dirname(path)
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) stop("Refusing to overwrite existing file: ", path)
  tmp <- tempfile(pattern = paste0(".", basename(path), ".tmp-"), tmpdir = dir_path)
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writeLines(as.character(text), tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) stop("Atomic rename failed for: ", path)
  invisible(path)
}

run_timed <- function(fun) {
  warnings <- character(0)
  start <- unname(proc.time()[["elapsed"]])
  value <- tryCatch(
    withCallingHandlers(
      fun(),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  elapsed <- unname(proc.time()[["elapsed"]]) - start
  if (inherits(value, "error")) {
    list(ok = FALSE, value = NULL, err = conditionMessage(value),
         warning = paste(unique(warnings), collapse = " | "), elapsed = elapsed)
  } else {
    list(ok = TRUE, value = value, err = "",
         warning = paste(unique(warnings), collapse = " | "), elapsed = elapsed)
  }
}

interp_Tmat_time_sensitivity <- function(transport_matrix, from_time, to_time) {
  transport_matrix <- as.matrix(transport_matrix)
  out <- matrix(NA_real_, nrow = length(to_time), ncol = ncol(transport_matrix))
  for (g_index in seq_len(ncol(transport_matrix))) {
    out[, g_index] <- approx(
      x = as.numeric(from_time), y = transport_matrix[, g_index],
      xout = as.numeric(to_time), rule = 2
    )$y
  }
  out
}

evaluate_transport_fit <- function(fit, truth, t_grid_out, xgrid) {
  estimated <- fit$T_hat
  if (length(fit$workGrid) != length(t_grid_out) ||
      max(abs(as.numeric(fit$workGrid) - t_grid_out)) > 1e-12) {
    estimated <- lapply(
      estimated,
      function(x) interp_Tmat_time_sensitivity(x, fit$workGrid, t_grid_out)
    )
  }

  n_subjects <- length(estimated)
  if (length(truth) != n_subjects) stop("Estimated/truth subject-count mismatch")
  imse_subject <- numeric(n_subjects)
  imse_sq_subject <- numeric(n_subjects)
  for (subject_index in seq_len(n_subjects)) {
    estimate_i <- as.matrix(estimated[[subject_index]])
    truth_i <- as.matrix(truth[[subject_index]])
    if (!all(dim(estimate_i) == dim(truth_i))) {
      stop("Transport dimension mismatch after time interpolation")
    }
    w1_t <- trapz_rows(xgrid, abs(estimate_i - truth_i))
    imse_subject[subject_index] <- trapz_vec(t_grid_out, w1_t)
    imse_sq_subject[subject_index] <- trapz_vec(t_grid_out, w1_t^2)
  }

  imse <- mean(imse_subject)
  imse_sq <- mean(imse_sq_subject)
  if (!is.finite(imse) || !is.finite(imse_sq)) {
    stop("Computed performance metric is non-finite")
  }
  list(imse = imse, imse_sq = imse_sq)
}

scalar_num_or_na <- function(x) {
  if (is.null(x) || length(x) < 1L) return(NA_real_)
  value <- suppressWarnings(as.numeric(x[1L]))
  if (length(value) != 1L || !is.finite(value)) NA_real_ else value
}

scalar_int_or_na <- function(x) {
  value <- scalar_num_or_na(x)
  if (is.na(value)) NA_integer_ else as.integer(value)
}

scalar_chr_or_na <- function(x) {
  if (is.null(x) || length(x) < 1L) return(NA_character_)
  as.character(x[1L])
}

extract_fit_info <- function(fit) {
  if (is.null(fit$fpca_fit) || is.null(fit$fpca_fit$optns)) {
    stop("fit_transport_fpca returned no fpca_fit/optns object")
  }
  fpca <- fit$fpca_fit
  optns <- fpca$optns
  diagnostics <- fit$diagnostics %||% list()
  internal_timings <- fpca$timings %||% numeric(0)

  list(
    fit_latent = if (is.logical(fit$latent) && length(fit$latent) == 1L) fit$latent else NA,
    actual_kernel = scalar_chr_or_na(optns$kernel),
    actual_bw_mu = scalar_num_or_na(fpca$bwMu),
    actual_bw_cov = scalar_num_or_na(fpca$bwCov),
    option_user_bw_mu = scalar_num_or_na(optns$userBwMu),
    option_user_bw_cov = scalar_num_or_na(optns$userBwCov),
    fpca_data_type = scalar_chr_or_na(optns$dataType),
    fpca_method_bw_mu = scalar_chr_or_na(optns$methodBwMu),
    fpca_method_bw_cov = scalar_chr_or_na(optns$methodBwCov),
    fpca_method_mu_cov_est = scalar_chr_or_na(optns$methodMuCovEst),
    fpca_method_xi = scalar_chr_or_na(optns$methodXi),
    fpca_use_binned_data = scalar_chr_or_na(optns$useBinnedData),
    fpca_use_binned_cov = if (is.null(optns$useBinnedCov)) NA else isTRUE(optns$useBinnedCov),
    fpca_n_reg_grid = scalar_int_or_na(optns$nRegGrid),
    fpca_select_k = scalar_int_or_na(fpca$selectK),
    fpca_fve = scalar_num_or_na(fpca$FVE),
    fpca_sigma2 = scalar_num_or_na(fpca$sigma2),
    fpca_rho = scalar_num_or_na(fpca$rho),
    fpca_internal_total_sec = scalar_num_or_na(internal_timings["total"]),
    fpca_internal_mu_sec = scalar_num_or_na(internal_timings["mu"]),
    fpca_internal_cov_sec = scalar_num_or_na(internal_timings["cov"]),
    fpca_internal_pace_sec = scalar_num_or_na(internal_timings["pace"]),
    diagnostic_n_obs_total = scalar_int_or_na(diagnostics$n_obs_total),
    diagnostic_kappa = scalar_num_or_na(diagnostics$kappa),
    diagnostic_u_cap = scalar_num_or_na(diagnostics$u_cap),
    diagnostic_max_abs_U_raw = scalar_num_or_na(diagnostics$max_abs_U_raw),
    diagnostic_frac_abs_U_raw_gt1 = scalar_num_or_na(diagnostics$frac_abs_U_raw_gt1),
    diagnostic_max_abs_U_used = scalar_num_or_na(diagnostics$max_abs_U_used),
    diagnostic_Z_cap = scalar_num_or_na(diagnostics$Z_cap),
    diagnostic_max_abs_Z_used = scalar_num_or_na(diagnostics$max_abs_Z_used)
  )
}

## These post-fit option-contract validators intentionally stop the whole run.
## They are called outside run_timed(): treating a contract mismatch as an
## ordinary failed fit would silently analyze a design other than the recorded
## kernel/bandwidth design.
validate_latent_fit <- function(info, latent, Z_cap, cell_label) {
  if (!identical(info$fit_latent, latent)) {
    stop("LATENT_OPTION_MISMATCH at ", cell_label)
  }
  if (latent) {
    if (!close_numeric(info$diagnostic_Z_cap, Z_cap) ||
        !is.finite(info$diagnostic_max_abs_Z_used) ||
        info$diagnostic_max_abs_Z_used < 0 ||
        info$diagnostic_max_abs_Z_used > Z_cap + 1e-10) {
      stop("LATENT_DIAGNOSTIC_MISMATCH at ", cell_label)
    }
  } else if (!is.na(info$diagnostic_Z_cap) ||
             !is.na(info$diagnostic_max_abs_Z_used)) {
    stop("NONLATENT_DIAGNOSTIC_MISMATCH at ", cell_label)
  }
  invisible(TRUE)
}

validate_default_fit <- function(info, time_range, cell_label, methodXi, latent, Z_cap) {
  validate_latent_fit(info, latent, Z_cap, cell_label)
  expected_bw_mu <- 0.05 * time_range
  expected_bw_cov <- 0.10 * time_range
  mismatches <- character(0)
  if (!identical(info$actual_kernel, "gauss")) {
    mismatches <- c(mismatches, paste0("kernel=", info$actual_kernel))
  }
  if (!identical(info$fpca_data_type, "Sparse")) {
    mismatches <- c(mismatches, paste0("dataType=", info$fpca_data_type))
  }
  if (!identical(info$fpca_method_mu_cov_est, "smooth")) {
    mismatches <- c(mismatches, paste0("methodMuCovEst=", info$fpca_method_mu_cov_est))
  }
  if (!identical(info$fpca_method_xi, methodXi)) {
    mismatches <- c(mismatches, paste0("methodXi=", info$fpca_method_xi))
  }
  if (!identical(info$fpca_method_bw_mu, "Default")) {
    mismatches <- c(mismatches, paste0("methodBwMu=", info$fpca_method_bw_mu))
  }
  if (!identical(info$fpca_method_bw_cov, "Default")) {
    mismatches <- c(mismatches, paste0("methodBwCov=", info$fpca_method_bw_cov))
  }
  if (!close_numeric(info$actual_bw_mu, expected_bw_mu)) {
    mismatches <- c(mismatches, sprintf("bwMu=%.17g expected %.17g",
                                        info$actual_bw_mu, expected_bw_mu))
  }
  if (!close_numeric(info$actual_bw_cov, expected_bw_cov)) {
    mismatches <- c(mismatches, sprintf("bwCov=%.17g expected %.17g",
                                        info$actual_bw_cov, expected_bw_cov))
  }
  if (length(mismatches) > 0L) {
    stop("FDAPACE_DEFAULT_MISMATCH at ", cell_label, ": ",
         paste(mismatches, collapse = "; "))
  }
  invisible(TRUE)
}

validate_explicit_fit <- function(info, kernel, requested_bw_mu, requested_bw_cov,
                                  cell_label, config_id, methodXi, latent, Z_cap) {
  validate_latent_fit(info, latent, Z_cap, paste(cell_label, config_id, sep = " / "))
  mismatches <- character(0)
  if (!identical(info$actual_kernel, kernel)) {
    mismatches <- c(mismatches, paste0("kernel=", info$actual_kernel,
                                       " requested ", kernel))
  }
  if (!close_numeric(info$actual_bw_mu, requested_bw_mu)) {
    mismatches <- c(mismatches, sprintf("bwMu=%.17g requested %.17g",
                                        info$actual_bw_mu, requested_bw_mu))
  }
  if (!close_numeric(info$actual_bw_cov, requested_bw_cov)) {
    mismatches <- c(mismatches, sprintf("bwCov=%.17g requested %.17g",
                                        info$actual_bw_cov, requested_bw_cov))
  }
  if (!close_numeric(info$option_user_bw_mu, requested_bw_mu) ||
      !close_numeric(info$option_user_bw_cov, requested_bw_cov)) {
    mismatches <- c(mismatches, "returned user bandwidth options differ from request")
  }
  if (!identical(info$fpca_data_type, "Sparse") ||
      !identical(info$fpca_method_mu_cov_est, "smooth") ||
      !identical(info$fpca_method_xi, methodXi)) {
    mismatches <- c(mismatches, "returned FPCA method options differ from fixed design")
  }
  if (length(mismatches) > 0L) {
    stop("FDAPACE_EXPLICIT_OPTION_MISMATCH at ", cell_label, " / ", config_id,
         ": ", paste(mismatches, collapse = "; "))
  }
  invisible(TRUE)
}

empty_baseline_norm_info <- function() {
  list(
    min = NA_real_, mean = NA_real_, max = NA_real_, sd = NA_real_,
    max_rel_dev = NA_real_
  )
}

empty_fit_info <- function() {
  list(
    fit_latent = NA,
    actual_kernel = NA_character_, actual_bw_mu = NA_real_, actual_bw_cov = NA_real_,
    option_user_bw_mu = NA_real_, option_user_bw_cov = NA_real_,
    fpca_data_type = NA_character_, fpca_method_bw_mu = NA_character_,
    fpca_method_bw_cov = NA_character_, fpca_method_mu_cov_est = NA_character_,
    fpca_method_xi = NA_character_, fpca_use_binned_data = NA_character_,
    fpca_use_binned_cov = NA, fpca_n_reg_grid = NA_integer_,
    fpca_select_k = NA_integer_, fpca_fve = NA_real_, fpca_sigma2 = NA_real_,
    fpca_rho = NA_real_, fpca_internal_total_sec = NA_real_,
    fpca_internal_mu_sec = NA_real_, fpca_internal_cov_sec = NA_real_,
    fpca_internal_pace_sec = NA_real_, diagnostic_n_obs_total = NA_integer_,
    diagnostic_kappa = NA_real_, diagnostic_u_cap = NA_real_,
    diagnostic_max_abs_U_raw = NA_real_, diagnostic_frac_abs_U_raw_gt1 = NA_real_,
    diagnostic_max_abs_U_used = NA_real_, diagnostic_Z_cap = NA_real_,
    diagnostic_max_abs_Z_used = NA_real_
  )
}

make_result_row <- function(cfg, tune, rep_id, n, N, m, seed,
                            base_bw_mu = NA_real_, base_bw_cov = NA_real_,
                            requested_bw_mu = NA_real_, requested_bw_cov = NA_real_,
                            fit_info = empty_fit_info(), imse = NA_real_, imse_sq = NA_real_,
                            ok = FALSE, err = "", stage = "unknown",
                            bandwidth_mode = "unavailable", user_bw_supplied = FALSE,
                            time_range = NA_real_, kappa_true = NA_real_,
                            baseline_norm_info = empty_baseline_norm_info(),
                            simulation_sec = NA_real_, truth_sec = NA_real_,
                            fit_sec = NA_real_, evaluation_sec = NA_real_,
                            data_warning = "", truth_warning = "", fit_warning = "") {
  core <- data.frame(
    run_id = cfg$run_id,
    setting = cfg$setting_id,
    rep = as.integer(rep_id),
    n = as.integer(n),
    N = as.integer(N),
    m = as.integer(m),
    seed = as.integer(seed),
    config_id = as.character(tune$config_id),
    kernel = as.character(tune$kernel),
    bw_multiplier = as.numeric(tune$bw_multiplier),
    is_default = isTRUE(tune$is_default),
    base_bw_mu = as.numeric(base_bw_mu),
    base_bw_cov = as.numeric(base_bw_cov),
    requested_bw_mu = as.numeric(requested_bw_mu),
    requested_bw_cov = as.numeric(requested_bw_cov),
    actual_bw_mu = as.numeric(fit_info$actual_bw_mu),
    actual_bw_cov = as.numeric(fit_info$actual_bw_cov),
    actual_kernel = as.character(fit_info$actual_kernel),
    imse = as.numeric(imse),
    imse_sq = as.numeric(imse_sq),
    ok = isTRUE(ok),
    err = as.character(err),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  extra <- data.frame(
    config_order = as.integer(tune$config_order),
    kernel_bw_scale = as.numeric(tune$kernel_bw_scale),
    mode = cfg$mode,
    setting_label = cfg$setting_label,
    status = if (isTRUE(ok)) "success" else "failure",
    stage = as.character(stage),
    bandwidth_mode = as.character(bandwidth_mode),
    user_bw_supplied = isTRUE(user_bw_supplied),
    latent = cfg$latent,
    fit_latent = as.logical(fit_info$fit_latent),
    methodXi = cfg$methodXi,
    methodMuCovEst = cfg$methodMuCovEst,
    eps_u = cfg$eps_u,
    u_cap = cfg$u_cap,
    Z_cap = cfg$Z_cap,
    decayR = cfg$decayR,
    beta_a_lower = cfg$beta_a[1L],
    beta_a_upper = cfg$beta_a[2L],
    beta_b_lower = cfg$beta_b[1L],
    beta_b_upper = cfg$beta_b[2L],
    t_design = cfg$t_design,
    t_grid_len = cfg$t_grid_len,
    xgrid_len = cfg$xgrid_len,
    generating_K = cfg$generating_K,
    seed_base = cfg$seed_base,
    pooled_time_range = as.numeric(time_range),
    kappa_true = as.numeric(kappa_true),
    baseline_norm_min = as.numeric(baseline_norm_info$min),
    baseline_norm_mean = as.numeric(baseline_norm_info$mean),
    baseline_norm_max = as.numeric(baseline_norm_info$max),
    baseline_norm_sd = as.numeric(baseline_norm_info$sd),
    baseline_norm_max_rel_dev = as.numeric(baseline_norm_info$max_rel_dev),
    option_user_bw_mu = as.numeric(fit_info$option_user_bw_mu),
    option_user_bw_cov = as.numeric(fit_info$option_user_bw_cov),
    fpca_data_type = as.character(fit_info$fpca_data_type),
    fpca_method_bw_mu = as.character(fit_info$fpca_method_bw_mu),
    fpca_method_bw_cov = as.character(fit_info$fpca_method_bw_cov),
    fpca_method_mu_cov_est = as.character(fit_info$fpca_method_mu_cov_est),
    fpca_method_xi = as.character(fit_info$fpca_method_xi),
    fpca_use_binned_data = as.character(fit_info$fpca_use_binned_data),
    fpca_use_binned_cov = as.logical(fit_info$fpca_use_binned_cov),
    fpca_n_reg_grid = as.integer(fit_info$fpca_n_reg_grid),
    fpca_select_k = as.integer(fit_info$fpca_select_k),
    fpca_fve = as.numeric(fit_info$fpca_fve),
    fpca_sigma2 = as.numeric(fit_info$fpca_sigma2),
    fpca_rho = as.numeric(fit_info$fpca_rho),
    diagnostic_n_obs_total = as.integer(fit_info$diagnostic_n_obs_total),
    diagnostic_kappa = as.numeric(fit_info$diagnostic_kappa),
    diagnostic_u_cap = as.numeric(fit_info$diagnostic_u_cap),
    diagnostic_max_abs_U_raw = as.numeric(fit_info$diagnostic_max_abs_U_raw),
    diagnostic_frac_abs_U_raw_gt1 = as.numeric(fit_info$diagnostic_frac_abs_U_raw_gt1),
    diagnostic_max_abs_U_used = as.numeric(fit_info$diagnostic_max_abs_U_used),
    diagnostic_Z_cap = as.numeric(fit_info$diagnostic_Z_cap),
    diagnostic_max_abs_Z_used = as.numeric(fit_info$diagnostic_max_abs_Z_used),
    simulation_sec = as.numeric(simulation_sec),
    truth_sec = as.numeric(truth_sec),
    fit_sec = as.numeric(fit_sec),
    evaluation_sec = as.numeric(evaluation_sec),
    fpca_internal_total_sec = as.numeric(fit_info$fpca_internal_total_sec),
    fpca_internal_mu_sec = as.numeric(fit_info$fpca_internal_mu_sec),
    fpca_internal_cov_sec = as.numeric(fit_info$fpca_internal_cov_sec),
    fpca_internal_pace_sec = as.numeric(fit_info$fpca_internal_pace_sec),
    data_warning = as.character(data_warning),
    truth_warning = as.character(truth_warning),
    fit_warning = as.character(fit_warning),
    fdapace_version = cfg$fdapace_version,
    dependency_md5 = cfg$dependency_md5,
    r_version = cfg$r_version,
    host = cfg$host,
    slurm_job_id = cfg$slurm_job_id,
    slurm_array_task_id = cfg$slurm_array_task_id,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cbind(core, extra)
}

run_one_rep <- function(rep_id, cfg, tuning_grid, rep_dir, mod_script) {
  rep_start <- unname(proc.time()[["elapsed"]])
  rows <- vector("list", length(cfg$n_vec) * length(cfg$N_vec) *
                   length(cfg$m_vec) * nrow(tuning_grid))
  row_index <- 0L
  default_index <- which(tuning_grid$is_default)
  default_tune <- tuning_grid[default_index, , drop = FALSE]

  append_row <- function(row) {
    row_index <<- row_index + 1L
    rows[[row_index]] <<- row
    invisible(NULL)
  }

  for (n in cfg$n_vec) {
    for (N in cfg$N_vec) {
      for (m in cfg$m_vec) {
        seed_num <- cfg$seed_base + rep_id * 1000000 + n * 1000 + N * 10 + m
        if (!is.finite(seed_num) || seed_num > .Machine$integer.max) {
          stop("Seed is outside R's integer range at rep=", rep_id,
               ", n=", n, ", N=", N, ", m=", m)
        }
        seed <- as.integer(seed_num)
        cell_label <- sprintf("rep=%d,n=%d,N=%d,m=%d,seed=%d", rep_id, n, N, m, seed)

        sim_result <- run_timed(function() {
          simulate_transport_process(
            n = n, N = N, m = m,
            decayR = cfg$decayR,
            xgrid = seq(0, 1, length.out = cfg$xgrid_len),
            K = cfg$generating_K,
            t_design = cfg$t_design,
            beta_a = cfg$beta_a,
            beta_b = cfg$beta_b,
            seed = seed,
            project = TRUE
          )
        })

        if (!sim_result$ok) {
          for (tune_index in seq_len(nrow(tuning_grid))) {
            tune <- tuning_grid[tune_index, , drop = FALSE]
            append_row(make_result_row(
              cfg, tune, rep_id, n, N, m, seed,
              ok = FALSE, err = paste0("Simulation failed: ", sim_result$err),
              stage = "simulation", simulation_sec = sim_result$elapsed,
              data_warning = sim_result$warning
            ))
          }
          next
        }

        sim <- sim_result$value
        baseline_norms <- transport_norm1_rows(sim$T0_true, sim$xgrid)
        if (length(baseline_norms) != n || any(!is.finite(baseline_norms)) ||
            !is.finite(sim$kappa_true) || sim$kappa_true <= 0) {
          stop("Invalid realized baseline-transport norms at ", cell_label)
        }
        baseline_norm_info <- list(
          min = min(baseline_norms),
          mean = mean(baseline_norms),
          max = max(baseline_norms),
          sd = if (length(baseline_norms) > 1L) stats::sd(baseline_norms) else 0,
          max_rel_dev = max(abs(baseline_norms - sim$kappa_true)) / abs(sim$kappa_true)
        )
        time_range <- diff(range(unlist(sim$t_list)))
        if (!is.finite(time_range) || time_range <= 0) {
          stop("Invalid pooled observation-time range at ", cell_label)
        }

        truth_result <- run_timed(function() {
          lapply(seq_len(n), function(subject_index) {
            true_transports_on_grid(
              Xi_row = sim$Xi[subject_index, ],
              T0_row = sim$T0_true[subject_index, ],
              t_grid = cfg$t_grid_out,
              xgrid = sim$xgrid,
              K = sim$K,
              g = sim$g,
              project = TRUE
            )
          })
        })

        if (!truth_result$ok) {
          for (tune_index in seq_len(nrow(tuning_grid))) {
            tune <- tuning_grid[tune_index, , drop = FALSE]
            append_row(make_result_row(
              cfg, tune, rep_id, n, N, m, seed,
              ok = FALSE, err = paste0("Truth evaluation failed: ", truth_result$err),
              stage = "truth", time_range = time_range,
              kappa_true = sim$kappa_true,
              baseline_norm_info = baseline_norm_info,
              simulation_sec = sim_result$elapsed, truth_sec = truth_result$elapsed,
              data_warning = sim_result$warning, truth_warning = truth_result$warning
            ))
          }
          rm(sim)
          gc(verbose = FALSE)
          next
        }
        truth <- truth_result$value

        fit_with_options <- function(fpca_optns) {
          fit_transport_fpca(
            T_list = sim$T_hat,
            t_list = sim$t_list,
            xgrid = sim$xgrid,
            kappa = sim$kappa_true,
            latent = cfg$latent,
            methodXi = cfg$methodXi,
            fpca_optns = fpca_optns,
            eps_u = cfg$eps_u,
            u_cap = cfg$u_cap,
            Z_cap = cfg$Z_cap,
            project = TRUE,
            verbose = FALSE
          )
        }

        ## The reference call deliberately contains no kernel or user bandwidth.
        default_fit_result <- run_timed(function() {
          fit_with_options(list(methodMuCovEst = cfg$methodMuCovEst))
        })

        if (!default_fit_result$ok) {
          for (tune_index in seq_len(nrow(tuning_grid))) {
            tune <- tuning_grid[tune_index, , drop = FALSE]
            error_prefix <- if (isTRUE(tune$is_default)) {
              "Implicit default fit failed: "
            } else {
              "Reference bandwidth unavailable because implicit default fit failed: "
            }
            append_row(make_result_row(
              cfg, tune, rep_id, n, N, m, seed,
              ok = FALSE, err = paste0(error_prefix, default_fit_result$err),
              stage = if (isTRUE(tune$is_default)) "fit_default" else "reference_unavailable",
              time_range = time_range, kappa_true = sim$kappa_true,
              baseline_norm_info = baseline_norm_info,
              simulation_sec = sim_result$elapsed, truth_sec = truth_result$elapsed,
              fit_sec = if (isTRUE(tune$is_default)) default_fit_result$elapsed else NA_real_,
              data_warning = sim_result$warning, truth_warning = truth_result$warning,
              fit_warning = if (isTRUE(tune$is_default)) default_fit_result$warning else ""
            ))
          }
          rm(sim, truth)
          gc(verbose = FALSE)
          next
        }

        default_fit <- default_fit_result$value
        default_info <- extract_fit_info(default_fit)
        validate_default_fit(default_info, time_range, cell_label,
                             cfg$methodXi, cfg$latent, cfg$Z_cap)
        base_bw_mu <- default_info$actual_bw_mu
        base_bw_cov <- default_info$actual_bw_cov

        default_eval_result <- run_timed(function() {
          evaluate_transport_fit(default_fit, truth, cfg$t_grid_out, sim$xgrid)
        })
        append_row(make_result_row(
          cfg, default_tune, rep_id, n, N, m, seed,
          base_bw_mu = base_bw_mu, base_bw_cov = base_bw_cov,
          requested_bw_mu = base_bw_mu, requested_bw_cov = base_bw_cov,
          fit_info = default_info,
          imse = if (default_eval_result$ok) default_eval_result$value$imse else NA_real_,
          imse_sq = if (default_eval_result$ok) default_eval_result$value$imse_sq else NA_real_,
          ok = default_eval_result$ok,
          err = if (default_eval_result$ok) "" else paste0("Default-fit evaluation failed: ", default_eval_result$err),
          stage = if (default_eval_result$ok) "complete" else "evaluation",
          bandwidth_mode = "implicit_default", user_bw_supplied = FALSE,
          time_range = time_range, kappa_true = sim$kappa_true,
          baseline_norm_info = baseline_norm_info,
          simulation_sec = sim_result$elapsed, truth_sec = truth_result$elapsed,
          fit_sec = default_fit_result$elapsed, evaluation_sec = default_eval_result$elapsed,
          data_warning = sim_result$warning, truth_warning = truth_result$warning,
          fit_warning = default_fit_result$warning
        ))
        default_fit_result$value <- NULL
        rm(default_fit)

        for (tune_index in setdiff(seq_len(nrow(tuning_grid)), default_index)) {
          tune <- tuning_grid[tune_index, , drop = FALSE]
          ## Match kernels by the leading local-linear second-moment bias scale before applying
          ## the prespecified perturbation around the Gaussian reference.
          requested_bw_mu <- as.numeric(tune$kernel_bw_scale) *
            as.numeric(tune$bw_multiplier) * base_bw_mu
          requested_bw_cov <- as.numeric(tune$kernel_bw_scale) *
            as.numeric(tune$bw_multiplier) * base_bw_cov
          explicit_optns <- list(
            methodMuCovEst = cfg$methodMuCovEst,
            kernel = as.character(tune$kernel),
            userBwMu = requested_bw_mu,
            userBwCov = requested_bw_cov
          )
          fit_result <- run_timed(function() fit_with_options(explicit_optns))

          if (!fit_result$ok) {
            append_row(make_result_row(
              cfg, tune, rep_id, n, N, m, seed,
              base_bw_mu = base_bw_mu, base_bw_cov = base_bw_cov,
              requested_bw_mu = requested_bw_mu, requested_bw_cov = requested_bw_cov,
              ok = FALSE, err = paste0("Explicit sensitivity fit failed: ", fit_result$err),
              stage = "fit", bandwidth_mode = "explicit_second_moment_scaled",
              user_bw_supplied = TRUE, time_range = time_range,
              kappa_true = sim$kappa_true,
              baseline_norm_info = baseline_norm_info,
              simulation_sec = sim_result$elapsed, truth_sec = truth_result$elapsed,
              fit_sec = fit_result$elapsed, data_warning = sim_result$warning,
              truth_warning = truth_result$warning, fit_warning = fit_result$warning
            ))
            next
          }

          fit <- fit_result$value
          fit_info <- extract_fit_info(fit)
          validate_explicit_fit(
            fit_info, as.character(tune$kernel), requested_bw_mu, requested_bw_cov,
            cell_label, as.character(tune$config_id), cfg$methodXi, cfg$latent, cfg$Z_cap
          )
          eval_result <- run_timed(function() {
            evaluate_transport_fit(fit, truth, cfg$t_grid_out, sim$xgrid)
          })
          append_row(make_result_row(
            cfg, tune, rep_id, n, N, m, seed,
            base_bw_mu = base_bw_mu, base_bw_cov = base_bw_cov,
            requested_bw_mu = requested_bw_mu, requested_bw_cov = requested_bw_cov,
            fit_info = fit_info,
            imse = if (eval_result$ok) eval_result$value$imse else NA_real_,
            imse_sq = if (eval_result$ok) eval_result$value$imse_sq else NA_real_,
            ok = eval_result$ok,
            err = if (eval_result$ok) "" else paste0("Sensitivity-fit evaluation failed: ", eval_result$err),
            stage = if (eval_result$ok) "complete" else "evaluation",
            bandwidth_mode = "explicit_second_moment_scaled", user_bw_supplied = TRUE,
            time_range = time_range, kappa_true = sim$kappa_true,
            baseline_norm_info = baseline_norm_info,
            simulation_sec = sim_result$elapsed, truth_sec = truth_result$elapsed,
            fit_sec = fit_result$elapsed, evaluation_sec = eval_result$elapsed,
            data_warning = sim_result$warning, truth_warning = truth_result$warning,
            fit_warning = fit_result$warning
          ))
          fit_result$value <- NULL
          rm(fit)
          gc(verbose = FALSE)
        }

        rm(sim, truth, default_fit_result, default_info)
        gc(verbose = FALSE)
      }
    }
  }

  rows <- rows[seq_len(row_index)]
  rep_df <- do.call(rbind, rows)
  expected_rep_rows <- length(cfg$n_vec) * length(cfg$N_vec) *
    length(cfg$m_vec) * nrow(tuning_grid)
  if (nrow(rep_df) != expected_rep_rows) {
    stop("Internal row-count mismatch for rep ", rep_id, ": got ", nrow(rep_df),
         ", expected ", expected_rep_rows)
  }
  rep_df$rep_elapsed_sec <- unname(proc.time()[["elapsed"]]) - rep_start
  rep_path <- file.path(rep_dir, sprintf("rep_%03d.rds", rep_id))
  atomic_save_rds(rep_df, rep_path)
  rep_df
}

## -------------------------------------------------------------------------
## Parse and validate run configuration.
## -------------------------------------------------------------------------

script_dir <- get_script_dir()
kv <- parse_kv_args(commandArgs(trailingOnly = TRUE))

allowed_keys <- c(
  "mode", "run_id", "outdir", "reps", "seed_base", "n_vec", "N_vec", "m_vec",
  "cores", "kernels", "bw_multipliers", "setting", "latent", "methodXi",
  "decayR", "beta_a", "beta_b", "t_design", "eps_u", "u_cap", "Z_cap",
  "t_grid_len", "mod_script", "check_only"
)
unknown_keys <- setdiff(names(kv), allowed_keys)
if (length(unknown_keys) > 0L) {
  stop("Unknown command-line key(s): ", paste(unknown_keys, collapse = ", "))
}

planned_n <- c(20L, 50L, 100L, 200L, 400L)
planned_N <- c(3L, 5L, 10L, 20L)
planned_m <- 200L
planned_kernels <- c("gauss", "epan", "quar")
planned_multipliers <- c(0.50, 0.75, 1.00, 1.25, 1.50)
kernel_second_moments <- c(gauss = 1, epan = 1 / 5, quar = 1 / 7)
kernel_bw_scales <- sqrt(kernel_second_moments[["gauss"]] / kernel_second_moments)
accepted_fdapace_kernels <- c("rect", "gauss", "epan", "gausvar", "quar")
expected_fdapace_version <- "0.6.0"
expected_r_version <- "4.5.0"
expected_dependency_md5 <- "4beac3d60a2b3e8850037ea7993abc05"

mode <- tolower(trimws(kv$mode %||% "full"))
if (!(mode %in% c("full", "smoke"))) stop("mode must be full or smoke")
check_only <- parse_bool(kv$check_only, FALSE, "check_only")

setting_arg <- tolower(gsub("[[:space:]_-]+", "", kv$setting %||% "1"))
if (!(setting_arg %in% c("1", "setting1"))) {
  stop("This prespecified driver implements only manuscript Setting 1")
}
setting_id <- 1L
setting_label <- "Setting 1"

reps <- parse_int_scalar(kv$reps, if (mode == "full") 200L else 1L,
                         "reps", lower = 1L)
seed_base <- parse_int_scalar(kv$seed_base, 2000L, "seed_base", lower = 1L)
n_vec <- parse_int_vec(kv$n_vec, if (mode == "full") planned_n else 20L, "n_vec")
N_vec <- parse_int_vec(kv$N_vec, if (mode == "full") planned_N else 3L, "N_vec")
m_vec <- parse_int_vec(kv$m_vec, planned_m, "m_vec")
kernels <- parse_char_vec(kv$kernels, planned_kernels, "kernels")
bw_multipliers <- parse_num_vec(kv$bw_multipliers, planned_multipliers,
                                "bw_multipliers")

if (!all(n_vec %in% planned_n) || !all(N_vec %in% planned_N) ||
    !all(m_vec %in% planned_m)) {
  stop("Smoke subsets must be drawn from the prespecified full n/N/m grid")
}
if (!all(kernels %in% planned_kernels)) {
  stop("kernels must be a subset of: ", paste(planned_kernels, collapse = ","))
}
if (!all(vapply(bw_multipliers, function(x) {
  any(vapply(planned_multipliers, function(y) close_numeric(x, y), logical(1L)))
}, logical(1L)))) {
  stop("bw_multipliers must be a subset of: ", paste(planned_multipliers, collapse = ","))
}
if (!("gauss" %in% kernels) || !any(vapply(bw_multipliers, function(x) {
  close_numeric(x, 1)
}, logical(1L)))) {
  stop("Every run, including smoke tests, must contain the gauss x 1.00 reference")
}

if (mode == "full") {
  if (reps != 200L || !identical(n_vec, planned_n) || !identical(N_vec, planned_N) ||
      !identical(m_vec, planned_m) || !identical(kernels, planned_kernels) ||
      !same_numeric_values(bw_multipliers, planned_multipliers)) {
    stop("mode=full requires exactly the prespecified 200-replication, 20-cell, 15-configuration design")
  }
}

latent <- parse_bool(kv$latent, FALSE, "latent")
methodXi <- trimws(kv$methodXi %||% "IN")
if (!(methodXi %in% c("CE", "IN"))) stop("methodXi must be CE or IN")
if (latent && methodXi != "CE") stop("latent=TRUE requires methodXi=CE in this study")
decayR <- parse_num_scalar(kv$decayR, 1, "decayR", lower = 0, lower_open = TRUE)
beta_a <- parse_num_vec(kv$beta_a, c(0.6, 3), "beta_a")
beta_b <- parse_num_vec(kv$beta_b, c(0.6, 3), "beta_b")
t_design <- tolower(trimws(kv$t_design %||% "random"))
eps_u <- parse_num_scalar(kv$eps_u, 0.05, "eps_u", lower = 0, upper = 1,
                          lower_open = TRUE, upper_open = TRUE)
u_cap <- parse_num_scalar(kv$u_cap, 0.95, "u_cap", lower = 0, upper = 1,
                          lower_open = TRUE, upper_open = TRUE)
Z_cap <- parse_num_scalar(kv$Z_cap, 8, "Z_cap", lower = 0, lower_open = TRUE)
t_grid_len <- parse_int_scalar(kv$t_grid_len, 51L, "t_grid_len", lower = 2L)

if (!close_numeric(decayR, 1) ||
    !same_numeric_values(beta_a, c(0.6, 3)) ||
    !same_numeric_values(beta_b, c(0.6, 3)) || t_design != "random" ||
    !close_numeric(eps_u, 0.05) || !close_numeric(u_cap, 0.95) ||
    !close_numeric(Z_cap, 8) || t_grid_len != 51L || seed_base != 2000L) {
  stop("Scientific fixed factors differ from the prespecified Setting-1 design")
}

slurm_job_id <- Sys.getenv("SLURM_JOB_ID", unset = "")
slurm_array_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
if (!check_only && !nzchar(slurm_job_id)) {
  stop("Simulation execution is permitted only inside a Slurm job; use check_only=TRUE for local preflight")
}

allocated_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")))
if (!is.finite(allocated_cores) || allocated_cores < 1L) allocated_cores <- NA_integer_
default_cores <- if (is.na(allocated_cores)) 1L else allocated_cores
cores_requested <- parse_int_scalar(kv$cores, default_cores, "cores", lower = 1L)
if (!is.na(allocated_cores) && cores_requested > allocated_cores) {
  stop("cores exceeds SLURM_CPUS_PER_TASK allocation")
}
cores_used <- min(cores_requested, reps)

run_id_default <- if (nzchar(slurm_job_id)) {
  paste0(mode, if (latent) "_latent_CE" else if (methodXi == "IN") "_IN" else "",
         "_", slurm_job_id,
         if (nzchar(slurm_array_task_id)) paste0("_", slurm_array_task_id) else "")
} else {
  paste0("check_", format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid())
}
run_id <- trimws(kv$run_id %||% run_id_default)
if (!grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", run_id)) {
  stop("run_id may contain only letters, digits, period, underscore, and hyphen")
}

mod_script_arg <- kv$mod_script %||%
  file.path(script_dir, "..", "parallel_run", "transport_fpca_modular.R")
if (!file.exists(mod_script_arg) && !grepl("^/", mod_script_arg)) {
  mod_script_arg <- file.path(script_dir, mod_script_arg)
}
if (!file.exists(mod_script_arg)) stop("Cannot find modular dependency: ", mod_script_arg)
mod_script <- normalizePath(mod_script_arg, mustWork = TRUE)
dependency_md5 <- unname(tools::md5sum(mod_script))
if (!identical(dependency_md5, expected_dependency_md5)) {
  stop("Dependency checksum mismatch for ", mod_script, "; expected ",
       expected_dependency_md5, ", got ", dependency_md5)
}

r_version <- as.character(getRversion())
if (!identical(r_version, expected_r_version)) {
  stop("R version mismatch: expected ", expected_r_version, ", got ", r_version)
}

if (!requireNamespace("fdapace", quietly = TRUE)) {
  stop("Package 'fdapace' is required")
}
fdapace_version <- as.character(utils::packageVersion("fdapace"))
if (!identical(fdapace_version, expected_fdapace_version)) {
  stop("fdapace version mismatch: expected ", expected_fdapace_version,
       ", got ", fdapace_version)
}

module_probe <- new.env(parent = globalenv())
sys.source(mod_script, envir = module_probe)
required_module_functions <- c(
  "fit_transport_fpca", "simulate_transport_process", "true_transports_on_grid",
  "trapz_vec", "trapz_rows", "transport_norm1_rows"
)
missing_module_functions <- required_module_functions[
  !vapply(required_module_functions, exists, logical(1L),
          envir = module_probe, inherits = FALSE)
]
if (length(missing_module_functions) > 0L) {
  stop("Modular dependency is missing required function(s): ",
       paste(missing_module_functions, collapse = ", "))
}
required_fit_args <- c("T_list", "t_list", "xgrid", "latent", "methodXi", "fpca_optns")
missing_fit_args <- setdiff(
  required_fit_args,
  names(formals(get("fit_transport_fpca", envir = module_probe, inherits = FALSE)))
)
if (length(missing_fit_args) > 0L) {
  stop("fit_transport_fpca lacks required argument(s): ", paste(missing_fit_args, collapse = ", "))
}
rm(module_probe)

## Probe the actual fdapace 0.6.0 option resolver without fitting or simulating.
set_options <- getFromNamespace("SetOptions", "fdapace")
probe_t <- list(
  c(0.11, 0.43, 0.81), c(0.19, 0.52, 0.93), c(0.14, 0.48, 0.86),
  c(0.23, 0.57, 0.89), c(0.16, 0.46, 0.84)
)
probe_y <- list(
  c(-0.4, 0.1, 0.5), c(-0.2, 0.2, 0.6), c(-0.3, 0.0, 0.4),
  c(-0.1, 0.3, 0.7), c(-0.5, -0.1, 0.3)
)
implicit_probe <- set_options(
  probe_y, probe_t, list(methodMuCovEst = "smooth", methodXi = methodXi)
)
probe_range <- diff(range(unlist(probe_t)))
if (!identical(implicit_probe$kernel, "gauss") ||
    !identical(implicit_probe$dataType, "Sparse") ||
    !identical(implicit_probe$methodXi, methodXi) ||
    !identical(implicit_probe$methodBwMu, "Default") ||
    !identical(implicit_probe$methodBwCov, "Default") ||
    !close_numeric(implicit_probe$userBwMu, 0.05 * probe_range) ||
    !close_numeric(implicit_probe$userBwCov, 0.10 * probe_range)) {
  stop("FDAPACE_DEFAULT_MISMATCH in package preflight")
}
explicit_probe <- set_options(
  probe_y, probe_t,
  list(
    methodMuCovEst = "smooth",
    methodXi = methodXi,
    kernel = implicit_probe$kernel,
    userBwMu = implicit_probe$userBwMu,
    userBwCov = implicit_probe$userBwCov
  )
)
critical_probe_fields <- c(
  "kernel", "userBwMu", "userBwCov", "methodBwMu", "methodBwCov",
  "methodMuCovEst", "dataType", "methodXi"
)
if (!isTRUE(all.equal(
  implicit_probe[critical_probe_fields], explicit_probe[critical_probe_fields],
  tolerance = 0, check.attributes = TRUE
))) {
  stop("Explicit default kernel/user bandwidths do not reproduce implicit fdapace options")
}
explicit_default_option_equivalence <- TRUE
rm(set_options, probe_t, probe_y)

tuning_grid <- do.call(
  rbind,
  lapply(kernels, function(kernel_name) {
    data.frame(
      kernel = kernel_name,
      kernel_second_moment = unname(kernel_second_moments[[kernel_name]]),
      kernel_bw_scale = unname(kernel_bw_scales[[kernel_name]]),
      bw_multiplier = bw_multipliers,
      stringsAsFactors = FALSE
    )
  })
)
tuning_grid$config_id <- sprintf("%s_x%.2f", tuning_grid$kernel,
                                 tuning_grid$bw_multiplier)
tuning_grid$is_default <- tuning_grid$kernel == "gauss" &
  vapply(tuning_grid$bw_multiplier, function(x) close_numeric(x, 1), logical(1L))
tuning_grid$config_order <- seq_len(nrow(tuning_grid))
if (sum(tuning_grid$is_default) != 1L || anyDuplicated(tuning_grid$config_id)) {
  stop("Tuning grid does not contain exactly one unique gauss x 1.00 reference")
}
expected_tuning_scales <- unname(kernel_bw_scales[tuning_grid$kernel])
if (any(!is.finite(tuning_grid$kernel_bw_scale)) ||
    !isTRUE(all.equal(tuning_grid$kernel_bw_scale, expected_tuning_scales,
                      tolerance = 0, check.attributes = FALSE))) {
  stop("Internal kernel bandwidth-scale mapping mismatch")
}

t_grid_out <- seq(0, 1, length.out = t_grid_len)
procedure_dir <- if (latent) "estimator14_CE" else if (methodXi == "IN") "estimator13_IN" else "direct_U_CE"
outdir_arg <- kv$outdir %||% file.path(script_dir, "results",
  paste0(if (mode == "smoke") "smoke_" else "", procedure_dir))

cfg <- list(
  run_id = run_id,
  mode = mode,
  setting_id = setting_id,
  setting_label = setting_label,
  reps = reps,
  seed_base = seed_base,
  n_vec = n_vec,
  N_vec = N_vec,
  m_vec = m_vec,
  kernels = kernels,
  bw_multipliers = bw_multipliers,
  kernel_second_moments = kernel_second_moments,
  kernel_bw_scales = kernel_bw_scales,
  latent = latent,
  methodXi = methodXi,
  methodMuCovEst = "smooth",
  decayR = decayR,
  beta_a = beta_a,
  beta_b = beta_b,
  t_design = t_design,
  eps_u = eps_u,
  u_cap = u_cap,
  Z_cap = Z_cap,
  t_grid_len = t_grid_len,
  t_grid_out = t_grid_out,
  xgrid_len = 128L,
  generating_K = 50L,
  cores_requested = cores_requested,
  cores_used = cores_used,
  fdapace_version = fdapace_version,
  accepted_fdapace_kernels = accepted_fdapace_kernels,
  expected_r_version = expected_r_version,
  expected_fdapace_version = expected_fdapace_version,
  dependency_md5 = dependency_md5,
  expected_dependency_md5 = expected_dependency_md5,
  mod_script = mod_script,
  explicit_default_option_equivalence = explicit_default_option_equivalence,
  r_version = r_version,
  host = Sys.info()[["nodename"]] %||% "",
  slurm_job_id = slurm_job_id,
  slurm_array_task_id = slurm_array_task_id,
  planned_rows = reps * length(n_vec) * length(N_vec) * length(m_vec) *
    nrow(tuning_grid)
)

if (check_only) {
  cat("Sensitivity driver preflight passed.\n")
  cat("R version:", r_version, "\n")
  cat("fdapace version:", fdapace_version, "\n")
  cat("dependency md5:", dependency_md5, "\n")
  cat("procedure: latent=", latent, ", methodXi=", methodXi, "\n", sep = "")
  cat("implicit defaults: kernel=gauss, bwMu=0.05*time_range, bwCov=0.10*time_range\n")
  cat("standard kernel second moments:",
      paste(paste0(names(kernel_second_moments), "=", format(kernel_second_moments,
                                                               digits = 16)),
            collapse = ","), "\n")
  cat("second-moment-equivalent kernel bandwidth scales:",
      paste(paste0(names(kernel_bw_scales), "=", format(kernel_bw_scales, digits = 16)),
            collapse = ","), "\n")
  cat("requested h formula: base Gaussian h * kernel_bw_scale * bw_multiplier\n")
  cat("accepted fdapace kernels:", paste(accepted_fdapace_kernels, collapse = ","), "\n")
  cat("planned run mode:", mode, "; rows:", cfg$planned_rows, "\n")
  quit(save = "no", status = 0L)
}

if (dir.exists(outdir_arg)) {
  existing <- setdiff(list.files(outdir_arg, all.files = TRUE), c(".", ".."))
  allowed_batch_entries <- "batch_manifest.txt"
  unexpected_entries <- setdiff(existing, allowed_batch_entries)
  if (length(unexpected_entries) > 0L) {
    stop("Refusing to use output directory containing unexpected existing entries: ",
         outdir_arg, " [", paste(unexpected_entries, collapse = ", "), "]")
  }
} else if (!dir.create(outdir_arg, recursive = TRUE, showWarnings = FALSE)) {
  stop("Cannot create output directory: ", outdir_arg)
}
outdir <- normalizePath(outdir_arg, mustWork = TRUE)
rep_dir <- file.path(outdir, "replications")
if (!dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(rep_dir)) {
  stop("Cannot create replication directory: ", rep_dir)
}
cfg$outdir <- outdir

run_start <- Sys.time()
config_bundle <- list(
  config = cfg,
  tuning_grid = tuning_grid,
  package_default_probe = implicit_probe,
  explicit_default_probe = explicit_probe,
  created_at = run_start
)
atomic_save_rds(config_bundle, file.path(outdir, "run_config.rds"))

config_csv <- cbind(
  data.frame(
    run_id = cfg$run_id, mode = cfg$mode, setting = cfg$setting_id,
    setting_label = cfg$setting_label,
    reps = cfg$reps, n_vec = paste(cfg$n_vec, collapse = ","),
    N_vec = paste(cfg$N_vec, collapse = ","), m_vec = paste(cfg$m_vec, collapse = ","),
    seed_base = cfg$seed_base, latent = cfg$latent, methodXi = cfg$methodXi,
    methodMuCovEst = cfg$methodMuCovEst, decayR = cfg$decayR,
    beta_a = paste(cfg$beta_a, collapse = ","), beta_b = paste(cfg$beta_b, collapse = ","),
    t_design = cfg$t_design, eps_u = cfg$eps_u, u_cap = cfg$u_cap, Z_cap = cfg$Z_cap,
    t_grid_len = cfg$t_grid_len, xgrid_len = cfg$xgrid_len,
    generating_K = cfg$generating_K, cores_requested = cfg$cores_requested,
    cores_used = cfg$cores_used, fdapace_version = cfg$fdapace_version,
    dependency_md5 = cfg$dependency_md5,
    explicit_default_option_equivalence = cfg$explicit_default_option_equivalence,
    stringsAsFactors = FALSE
  )[rep(1L, nrow(tuning_grid)), , drop = FALSE],
  tuning_grid
)
atomic_write_csv(config_csv, file.path(outdir, "run_config.csv"))

cat("\n===== kernel/bandwidth sensitivity analysis =====\n")
cat("run_id:", cfg$run_id, "\n")
cat("mode:", cfg$mode, "\n")
cat("procedure: latent=", cfg$latent, ", methodXi=", cfg$methodXi, "\n", sep = "")
cat("outdir:", outdir, "\n")
cat("reps:", reps, "; cells per rep:",
    length(n_vec) * length(N_vec) * length(m_vec), "\n")
cat("configurations:", nrow(tuning_grid), "; planned rows:", cfg$planned_rows, "\n")
cat("cores:", cores_used, "of", cores_requested, "requested\n")
cat("fdapace:", fdapace_version, "; dependency md5:", dependency_md5, "\n")
cat("=================================================\n\n")

cluster <- parallel::makePSOCKcluster(cores_used, outfile = "")
on.exit({
  if (!is.null(cluster)) try(parallel::stopCluster(cluster), silent = TRUE)
}, add = TRUE)

parallel::clusterExport(
  cluster,
  varlist = c("mod_script", "expected_r_version", "expected_fdapace_version"),
  envir = environment()
)
invisible(parallel::clusterEvalQ(cluster, {
  suppressMessages(library(fdapace))
  if (as.character(getRversion()) != expected_r_version) {
    stop("Worker R version mismatch")
  }
  if (as.character(utils::packageVersion("fdapace")) != expected_fdapace_version) {
    stop("Worker fdapace version mismatch")
  }
  sys.source(mod_script, envir = .GlobalEnv)
  NULL
}))

worker_exports <- c(
  "%||%", "atomic_save_rds", "run_timed", "interp_Tmat_time_sensitivity",
  "evaluate_transport_fit", "scalar_num_or_na", "scalar_int_or_na",
  "scalar_chr_or_na", "extract_fit_info", "empty_fit_info", "close_numeric",
  "validate_latent_fit", "validate_default_fit", "validate_explicit_fit", "empty_baseline_norm_info",
  "make_result_row", "run_one_rep"
)
parallel::clusterExport(cluster, varlist = worker_exports, envir = environment())

rep_ids <- seq_len(reps)
simulation_clock <- proc.time()
rep_results <- parallel::parLapply(
  cluster, rep_ids, run_one_rep,
  cfg = cfg, tuning_grid = tuning_grid, rep_dir = rep_dir, mod_script = mod_script
)
run_time <- proc.time() - simulation_clock

parallel::stopCluster(cluster)
cluster <- NULL

raw_df <- do.call(rbind, rep_results)
rownames(raw_df) <- NULL
config_rank <- match(raw_df$config_id, tuning_grid$config_id)
raw_df <- raw_df[order(raw_df$rep, raw_df$n, raw_df$N, raw_df$m, config_rank), , drop = FALSE]
rownames(raw_df) <- NULL

core_columns <- c(
  "run_id", "setting", "rep", "n", "N", "m", "seed", "config_id", "kernel",
  "bw_multiplier", "is_default", "base_bw_mu", "base_bw_cov", "requested_bw_mu",
  "requested_bw_cov", "actual_bw_mu", "actual_bw_cov", "actual_kernel", "imse",
  "imse_sq", "ok", "err"
)
if (!identical(names(raw_df)[seq_along(core_columns)], core_columns)) {
  stop("Internal raw-result core schema mismatch")
}
if (nrow(raw_df) != cfg$planned_rows) {
  stop("Final raw row count is ", nrow(raw_df), "; expected ", cfg$planned_rows)
}
expected_row_scales <- unname(cfg$kernel_bw_scales[raw_df$kernel])
if (any(!is.finite(expected_row_scales)) ||
    any(!is.finite(raw_df$kernel_bw_scale)) ||
    !isTRUE(all.equal(raw_df$kernel_bw_scale, expected_row_scales,
                      tolerance = 0, check.attributes = FALSE))) {
  stop("Raw-result kernel bandwidth-scale mapping mismatch")
}

finite_base_mu <- is.finite(raw_df$base_bw_mu)
finite_base_cov <- is.finite(raw_df$base_bw_cov)
if (any(xor(finite_base_mu, finite_base_cov))) {
  stop("Raw-result Gaussian reference bandwidths are only partially available")
}
has_reference_bw <- finite_base_mu & finite_base_cov
if (any(!has_reference_bw &
        (is.finite(raw_df$requested_bw_mu) | is.finite(raw_df$requested_bw_cov)))) {
  stop("Raw-result requested bandwidth is present without Gaussian reference bandwidths")
}
if (any(has_reference_bw)) {
  expected_requested_mu <- raw_df$base_bw_mu[has_reference_bw] *
    raw_df$kernel_bw_scale[has_reference_bw] * raw_df$bw_multiplier[has_reference_bw]
  expected_requested_cov <- raw_df$base_bw_cov[has_reference_bw] *
    raw_df$kernel_bw_scale[has_reference_bw] * raw_df$bw_multiplier[has_reference_bw]
  observed_requested_mu <- raw_df$requested_bw_mu[has_reference_bw]
  observed_requested_cov <- raw_df$requested_bw_cov[has_reference_bw]
  mu_formula_ok <- is.finite(observed_requested_mu) &
    abs(observed_requested_mu - expected_requested_mu) <=
      1e-10 * pmax(1, abs(observed_requested_mu), abs(expected_requested_mu))
  cov_formula_ok <- is.finite(observed_requested_cov) &
    abs(observed_requested_cov - expected_requested_cov) <=
      1e-10 * pmax(1, abs(observed_requested_cov), abs(expected_requested_cov))
  if (!all(mu_formula_ok) || !all(cov_formula_ok)) {
    stop("Raw-result requested bandwidths violate base h * kernel_bw_scale * multiplier")
  }
}
key <- interaction(
  raw_df$rep, raw_df$n, raw_df$N, raw_df$m, raw_df$kernel,
  sprintf("%.8f", raw_df$bw_multiplier), drop = TRUE
)
if (anyDuplicated(key)) stop("Duplicate raw-result key detected")

cell_key <- interaction(raw_df$rep, raw_df$n, raw_df$N, raw_df$m, drop = TRUE)
cell_counts <- table(cell_key)
if (any(cell_counts != nrow(tuning_grid))) {
  stop("One or more generated-data cells do not have the planned number of configurations")
}
default_counts <- table(cell_key[raw_df$is_default])
if (length(default_counts) != length(cell_counts) || any(default_counts != 1L)) {
  stop("One or more generated-data cells lack exactly one default row")
}
cell_seed <- unique(raw_df[c("rep", "n", "N", "m", "seed")])
expected_cells <- reps * length(n_vec) * length(N_vec) * length(m_vec)
if (nrow(cell_seed) != expected_cells || anyDuplicated(cell_seed$seed)) {
  stop("Seed uniqueness/completeness validation failed")
}
baseline_columns <- c(
  "baseline_norm_min", "baseline_norm_mean", "baseline_norm_max",
  "baseline_norm_sd", "baseline_norm_max_rel_dev"
)
baseline_matrix <- as.matrix(raw_df[baseline_columns])
simulation_succeeded <- is.finite(raw_df$kappa_true)
finite_baseline_count <- rowSums(is.finite(baseline_matrix))
if (any(!simulation_succeeded & raw_df$stage != "simulation") ||
    any(simulation_succeeded & raw_df$stage == "simulation")) {
  stop("Simulation-stage labels are inconsistent with realized-data diagnostics")
}
if (any(simulation_succeeded & finite_baseline_count != length(baseline_columns)) ||
    any(!simulation_succeeded & finite_baseline_count != 0L)) {
  stop("Baseline-norm diagnostic availability is inconsistent with simulation status")
}
if (any(simulation_succeeded)) {
  sim_rows <- which(simulation_succeeded)
  norm_min <- raw_df$baseline_norm_min[sim_rows]
  norm_mean <- raw_df$baseline_norm_mean[sim_rows]
  norm_max <- raw_df$baseline_norm_max[sim_rows]
  norm_sd <- raw_df$baseline_norm_sd[sim_rows]
  norm_rel_dev <- raw_df$baseline_norm_max_rel_dev[sim_rows]
  kappa <- raw_df$kappa_true[sim_rows]
  expected_rel_dev <- pmax(abs(norm_min - kappa), abs(norm_max - kappa)) / abs(kappa)
  rel_dev_ok <- abs(norm_rel_dev - expected_rel_dev) <=
    1e-10 * pmax(1, abs(norm_rel_dev), abs(expected_rel_dev))
  if (any(kappa <= 0) || any(norm_min < 0) || any(norm_sd < 0) ||
      any(norm_rel_dev < 0) || any(norm_min > norm_mean) ||
      any(norm_mean > norm_max) || !all(rel_dev_ok)) {
    stop("Baseline-norm diagnostic structural validation failed")
  }
}
cell_diagnostics <- unique(raw_df[c(
  "rep", "n", "N", "m", "seed", "kappa_true", baseline_columns
)])
if (nrow(cell_diagnostics) != expected_cells) {
  stop("Baseline-norm diagnostics are not constant within generated-data cells")
}
if (any(raw_df$ok & (!is.finite(raw_df$imse) | !is.finite(raw_df$imse_sq)))) {
  stop("A successful row has a non-finite metric")
}
if (any(!raw_df$ok & !nzchar(raw_df$err))) {
  stop("A failed row has no error message")
}
rep_files <- list.files(rep_dir, pattern = "^rep_[0-9]+\\.rds$", full.names = TRUE)
if (length(rep_files) != reps) stop("Replication-file count mismatch")

raw_rds_path <- file.path(outdir, "raw_results.rds")
raw_csv_path <- file.path(outdir, "raw_results.csv")
atomic_save_rds(raw_df, raw_rds_path)
atomic_write_csv(raw_df, raw_csv_path)

session_lines <- capture.output(print(sessionInfo()))
atomic_write_lines(session_lines, file.path(outdir, "sessionInfo.txt"))

run_end <- Sys.time()
metadata <- list(
  run_id = cfg$run_id,
  status = if (all(raw_df$ok)) "complete" else "complete_with_failures",
  start_time = run_start,
  end_time = run_end,
  elapsed_wall_seconds = as.numeric(difftime(run_end, run_start, units = "secs")),
  parallel_proc_time = run_time,
  planned_rows = cfg$planned_rows,
  observed_rows = nrow(raw_df),
  successful_rows = sum(raw_df$ok),
  failed_rows = sum(!raw_df$ok),
  rep_files = length(rep_files),
  expected_cells = expected_cells,
  unique_seeds = length(unique(cell_seed$seed)),
  structural_validation_passed = TRUE,
  raw_results_rds_md5 = unname(tools::md5sum(raw_rds_path)),
  raw_results_csv_md5 = unname(tools::md5sum(raw_csv_path)),
  config = cfg,
  tuning_grid = tuning_grid
)
atomic_save_rds(metadata, file.path(outdir, "run_metadata.rds"))
metadata_text <- c(
  paste0("run_id: ", metadata$run_id),
  paste0("latent: ", cfg$latent),
  paste0("methodXi: ", cfg$methodXi),
  paste0("status: ", metadata$status),
  paste0("start_time: ", format(metadata$start_time, tz = "UTC", usetz = TRUE)),
  paste0("end_time: ", format(metadata$end_time, tz = "UTC", usetz = TRUE)),
  paste0("elapsed_wall_seconds: ", metadata$elapsed_wall_seconds),
  paste0("planned_rows: ", metadata$planned_rows),
  paste0("observed_rows: ", metadata$observed_rows),
  paste0("successful_rows: ", metadata$successful_rows),
  paste0("failed_rows: ", metadata$failed_rows),
  paste0("rep_files: ", metadata$rep_files),
  paste0("unique_seeds: ", metadata$unique_seeds),
  paste0("raw_results_rds_md5: ", metadata$raw_results_rds_md5),
  paste0("raw_results_csv_md5: ", metadata$raw_results_csv_md5)
)
atomic_write_lines(metadata_text, file.path(outdir, "run_metadata.txt"))

run_status <- data.frame(
  run_id = cfg$run_id,
  status = metadata$status,
  planned_rows = metadata$planned_rows,
  observed_rows = metadata$observed_rows,
  successful_rows = metadata$successful_rows,
  failed_rows = metadata$failed_rows,
  start_time_utc = format(run_start, tz = "UTC", usetz = TRUE),
  end_time_utc = format(run_end, tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
atomic_write_csv(run_status, file.path(outdir, "run_status.csv"))
atomic_write_lines(
  c(
    paste0("run_id=", cfg$run_id),
    paste0("status=", metadata$status),
    paste0("observed_rows=", metadata$observed_rows),
    paste0("failed_rows=", metadata$failed_rows)
  ),
  file.path(outdir, "RUN_COMPLETE")
)

cat("\n===== SENSITIVITY RUN COMPLETE =====\n")
cat("status:", metadata$status, "\n")
cat("rows:", metadata$observed_rows, "of", metadata$planned_rows, "\n")
cat("successful:", metadata$successful_rows, "; failed:", metadata$failed_rows, "\n")
cat("raw RDS:", raw_rds_path, "\n")
cat("raw CSV:", raw_csv_path, "\n")
cat("completion marker:", file.path(outdir, "RUN_COMPLETE"), "\n")
cat("====================================\n")
