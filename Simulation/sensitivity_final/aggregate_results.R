#!/usr/bin/env Rscript

## Deterministic validation and aggregation for the Setting-1 kernel/bandwidth
## sensitivity analysis.  This script uses base R only.  It never runs a
## simulation and refuses to reuse a nonempty aggregation directory.

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

get_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) > 0L) {
    script_path <- sub("^--file=", "", file_arg[1L])
    return(dirname(normalizePath(script_path, mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

parse_kv_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) stop("Every argument must have key=value form: ", arg)
    key <- sub("=.*$", "", arg)
    value <- sub("^[^=]*=", "", arg)
    if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", key) || !is.null(out[[key]])) {
      stop("Invalid or duplicated argument key: ", key)
    }
    out[[key]] <- value
  }
  out
}

split_csv <- function(x) strsplit(gsub("[[:space:]]+", "", x), ",", fixed = TRUE)[[1L]]
parse_int_vec <- function(x) as.integer(split_csv(x))
parse_num_vec <- function(x) as.numeric(split_csv(x))
parse_chr_vec <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1L]])

parse_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  value <- tolower(trimws(as.character(x)[1L]))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Cannot parse logical value: ", x)
}

normalize_logical <- function(x, name) {
  if (is.logical(x)) return(x)
  value <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(value))
  out[value %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[value %in% c("false", "f", "0", "no", "n")] <- FALSE
  if (anyNA(out)) stop("Column '", name, "' contains unparseable logical values.")
  out
}

close_vec <- function(x, y, tolerance = 1e-10) {
  x <- as.numeric(x)
  y <- rep(as.numeric(y), length.out = length(x))
  is.finite(x) & is.finite(y) & abs(x - y) <= tolerance * pmax(1, abs(x), abs(y))
}

same_numeric <- function(x, y, tolerance = 1e-12) {
  length(x) == length(y) && all(close_vec(x, y, tolerance))
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x) / sqrt(length(x))
}

numeric_key <- function(x) sprintf("%.8f", as.numeric(x))
dataset_key <- function(df) paste(df$run_id, df$setting, df$rep, df$n, df$N, df$m, sep = "|")
result_key <- function(df) paste(dataset_key(df), df$kernel, numeric_key(df$bw_multiplier), sep = "|")
kernel_order <- function(x) match(x, c("gauss", "epan", "quar"))

sort_results <- function(df, include_rep = FALSE) {
  if (!nrow(df)) return(df)
  args <- list(df$setting, kernel_order(df$kernel), df$bw_multiplier)
  for (name in c("n", "N", "m")) if (name %in% names(df)) args <- c(args, list(df[[name]]))
  if (include_rep && "rep" %in% names(df)) args <- c(args, list(df$rep))
  df[do.call(order, args), , drop = FALSE]
}

group_apply <- function(df, columns, fun) {
  if (!nrow(df)) return(data.frame())
  key <- do.call(paste, c(df[columns], sep = "\r"))
  pieces <- split(seq_len(nrow(df)), key)
  rows <- lapply(pieces, function(index) {
    cbind(df[index[1L], columns, drop = FALSE], fun(df[index, , drop = FALSE]))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

prepare_empty_dir <- function(path) {
  if (dir.exists(path)) {
    entries <- setdiff(list.files(path, all.files = TRUE), c(".", ".."))
    if (length(entries)) stop("Refusing to use nonempty aggregation directory: ", path)
  } else if (!dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create aggregation directory: ", path)
  }
  normalizePath(path, mustWork = TRUE)
}

atomic_write <- function(path, writer, allow_replace = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !allow_replace) stop("Refusing to overwrite existing file: ", path)
  temporary <- tempfile(pattern = paste0(".", basename(path), ".tmp-"), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  writer(temporary)
  if (!file.rename(temporary, path)) stop("Atomic rename failed for: ", path)
  invisible(path)
}

atomic_write_csv <- function(x, path) {
  atomic_write(path, function(tmp) utils::write.csv(x, tmp, row.names = FALSE, na = "NA"))
}

atomic_write_rds <- function(x, path) {
  atomic_write(path, function(tmp) saveRDS(x, tmp, version = 3))
}

atomic_write_lines <- function(x, path, allow_replace = FALSE) {
  atomic_write(path, function(tmp) writeLines(as.character(x), tmp, useBytes = TRUE), allow_replace)
}

kv <- parse_kv_args(commandArgs(trailingOnly = TRUE))
mode <- tolower(kv$mode %||% "full")
if (!mode %in% c("smoke", "full")) stop("mode must be smoke or full.")
strict <- parse_bool(kv$strict, TRUE)
write_table <- parse_bool(kv$write_table, mode == "full")

input_dir <- kv$input_dir %||% stop("Required argument: input_dir=PATH")
if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
input_dir <- normalizePath(input_dir, mustWork = TRUE)
outdir <- prepare_empty_dir(kv$outdir %||% file.path(input_dir, "aggregated"))
table_out <- kv$table_out %||% file.path(outdir, "sensitivity_table.tex")

expected_reps <- if (!is.null(kv$expected_reps)) as.integer(kv$expected_reps) else if (mode == "smoke") 2L else 200L
n_vec <- if (!is.null(kv$n_vec)) parse_int_vec(kv$n_vec) else if (mode == "smoke") 20L else c(20L, 50L, 100L, 200L, 400L)
N_vec <- if (!is.null(kv$N_vec)) parse_int_vec(kv$N_vec) else if (mode == "smoke") 3L else c(3L, 5L, 10L, 20L)
m_vec <- if (!is.null(kv$m_vec)) parse_int_vec(kv$m_vec) else 200L
kernels <- if (!is.null(kv$kernels)) parse_chr_vec(kv$kernels) else c("gauss", "epan", "quar")
bw_multipliers <- if (!is.null(kv$bw_multipliers)) parse_num_vec(kv$bw_multipliers) else c(0.50, 0.75, 1.00, 1.25, 1.50)

if (length(expected_reps) != 1L || !is.finite(expected_reps) || expected_reps < 1L ||
    any(!is.finite(n_vec)) || any(n_vec < 1L) || anyDuplicated(n_vec) ||
    any(!is.finite(N_vec)) || any(N_vec < 2L) || anyDuplicated(N_vec) ||
    any(!is.finite(m_vec)) || any(m_vec < 1L) || anyDuplicated(m_vec) ||
    !identical(kernels, c("gauss", "epan", "quar")) ||
    !same_numeric(bw_multipliers, c(0.50, 0.75, 1.00, 1.25, 1.50))) {
  stop("Requested aggregation design is not the fixed sensitivity design.")
}
if (mode == "full" && (!identical(as.integer(expected_reps), 200L) ||
    !identical(as.integer(n_vec), c(20L, 50L, 100L, 200L, 400L)) ||
    !identical(as.integer(N_vec), c(3L, 5L, 10L, 20L)) ||
    !identical(as.integer(m_vec), 200L))) {
  stop("Full mode requires exactly 200 repetitions and all 20 prespecified (n,N) cells.")
}
if (mode == "smoke" && (!identical(as.integer(n_vec), 20L) ||
    !identical(as.integer(N_vec), 3L) || !identical(as.integer(m_vec), 200L))) {
  stop("Smoke mode is restricted to the prespecified hardest cell n=20, N=3, m=200.")
}

checks <- list()
numerical_status <- "NOT_ASSESSED"
add_check <- function(name, passed, details) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, passed = isTRUE(passed), details = as.character(details),
    stringsAsFactors = FALSE
  )
}

write_validation <- function() {
  check_df <- if (length(checks)) do.call(rbind, checks) else data.frame(
    check = character(), passed = logical(), details = character(), stringsAsFactors = FALSE
  )
  rownames(check_df) <- NULL
  atomic_write_csv(check_df, file.path(outdir, "validation_checks.csv"))
  overall <- nrow(check_df) > 0L && all(check_df$passed)
  lines <- c(
    "Sensitivity analysis validation",
    paste0("Mode: ", mode),
    paste0("Input: ", input_dir),
    paste0("Structural validation: ", if (overall) "PASS" else "FAIL"),
    paste0("Numerical completion: ", numerical_status),
    "",
    if (nrow(check_df)) vapply(seq_len(nrow(check_df)), function(i) {
      sprintf("[%s] %s: %s", if (check_df$passed[i]) "PASS" else "FAIL",
              check_df$check[i], check_df$details[i])
    }, character(1L)) else character()
  )
  atomic_write_lines(lines, file.path(outdir, "validation_report.txt"))
  overall
}

## Read and reconcile the combined result with every atomic replication file.
expected_rep_names <- sprintf("rep_%03d.rds", seq_len(expected_reps))
rep_dir <- file.path(input_dir, "replications")
observed_rep_names <- if (dir.exists(rep_dir)) {
  sort(list.files(rep_dir, pattern = "^rep_[0-9]{3}\\.rds$", full.names = FALSE))
} else character()
add_check("replication_directory", dir.exists(rep_dir), rep_dir)
add_check("replication_file_names", identical(observed_rep_names, expected_rep_names),
          sprintf("Observed %d exact names; expected %d.", length(observed_rep_names), expected_reps))

rep_paths <- file.path(rep_dir, expected_rep_names)
rep_objects <- NULL
rep_read_error <- ""
if (all(file.exists(rep_paths))) {
  rep_objects <- tryCatch(lapply(rep_paths, readRDS), error = function(e) {
    rep_read_error <<- conditionMessage(e)
    NULL
  })
}
rep_read_ok <- !is.null(rep_objects) && length(rep_objects) == expected_reps &&
  all(vapply(rep_objects, is.data.frame, logical(1L)))
add_check("replication_files_readable", rep_read_ok,
          if (rep_read_ok) "All expected replication files contain data frames." else rep_read_error %||% "Replication files are missing or invalid.")

minimal_rep_columns <- c("rep", "n", "N", "m", "kernel", "bw_multiplier")
expected_rows_per_rep <- length(n_vec) * length(N_vec) * length(m_vec) * 15L
rep_identity_ok <- rep_read_ok && all(vapply(seq_len(expected_reps), function(i) {
  d <- rep_objects[[i]]
  all(minimal_rep_columns %in% names(d)) && nrow(d) == expected_rows_per_rep &&
    all(suppressWarnings(as.integer(d$rep)) == i) &&
    !anyDuplicated(paste(d$n, d$N, d$m, d$kernel, numeric_key(d$bw_multiplier), sep = "|"))
}, logical(1L)))
add_check("replication_file_identity", rep_identity_ok,
          sprintf("Each rep_### file must contain exactly %d unique rows for its filename replication.", expected_rows_per_rep))

rep_raw <- NULL
if (rep_read_ok) {
  rep_raw <- do.call(rbind, rep_objects)
  rownames(rep_raw) <- NULL
}

combined_path <- file.path(input_dir, "raw_results.rds")
combined_present <- file.exists(combined_path)
combined <- if (combined_present) tryCatch(readRDS(combined_path), error = identity) else NULL
combined_ok <- is.data.frame(combined)
add_check("combined_raw_readable_or_absent", combined_ok || (!combined_present && rep_read_ok),
          if (combined_ok) combined_path else if (!combined_present && rep_read_ok)
            "raw_results.rds is absent; reconstruction from complete replication files is allowed."
          else "raw_results.rds exists but is unreadable.")

if (combined_ok) rownames(combined) <- NULL
canonicalize <- function(d) {
  if (!is.data.frame(d) || !all(c("run_id", "setting", "rep", "n", "N", "m", "kernel", "bw_multiplier") %in% names(d))) return(d)
  d <- d[order(as.character(d$run_id), as.numeric(d$setting), as.numeric(d$rep),
               as.numeric(d$n), as.numeric(d$N), as.numeric(d$m),
               kernel_order(as.character(d$kernel)), as.numeric(d$bw_multiplier)), , drop = FALSE]
  rownames(d) <- NULL
  d
}
combined_match <- (!combined_present && rep_read_ok) || (combined_ok && rep_read_ok &&
  isTRUE(all.equal(canonicalize(combined), canonicalize(rep_raw),
                   tolerance = 0, check.attributes = TRUE)))
add_check("combined_matches_replications", combined_match,
          if (combined_match && combined_present) "Combined and replication results agree exactly after canonical key ordering."
          else if (combined_match) "Complete replication files are the canonical fallback input."
          else "Combined and replication results differ or cannot be compared.")

raw <- if (combined_ok) combined else rep_raw
if (!is.data.frame(raw)) {
  write_validation()
  stop("No readable raw result is available; see validation_report.txt.", call. = FALSE)
}
rownames(raw) <- NULL

source_files <- c(if (combined_ok) combined_path else character(), rep_paths[file.exists(rep_paths)])
atomic_write_lines(normalizePath(source_files, mustWork = TRUE), file.path(outdir, "input_files.txt"))

config_path <- file.path(input_dir, "run_config.rds")
config_bundle <- if (file.exists(config_path)) tryCatch(readRDS(config_path), error = identity) else NULL
config_readable <- is.list(config_bundle) && is.list(config_bundle$config) && is.data.frame(config_bundle$tuning_grid)
add_check("run_config_readable", config_readable, config_path)

metadata_path <- file.path(input_dir, "run_metadata.rds")
metadata <- if (file.exists(metadata_path)) tryCatch(readRDS(metadata_path), error = identity) else NULL
metadata_readable <- is.list(metadata)
add_check("run_metadata_readable", metadata_readable, metadata_path)
add_check("completion_marker", file.exists(file.path(input_dir, "RUN_COMPLETE")),
          file.path(input_dir, "RUN_COMPLETE"))

required_columns <- c(
  "run_id", "setting", "rep", "n", "N", "m", "seed", "config_id", "kernel",
  "bw_multiplier", "is_default", "base_bw_mu", "base_bw_cov", "requested_bw_mu",
  "requested_bw_cov", "actual_bw_mu", "actual_bw_cov", "actual_kernel", "imse",
  "imse_sq", "ok", "err", "config_order", "kernel_bw_scale", "mode",
  "setting_label", "status", "stage", "bandwidth_mode", "user_bw_supplied",
  "latent", "methodXi", "methodMuCovEst", "eps_u", "u_cap", "Z_cap", "decayR",
  "beta_a_lower", "beta_a_upper", "beta_b_lower", "beta_b_upper", "t_design",
  "t_grid_len", "xgrid_len", "generating_K", "seed_base", "pooled_time_range",
  "kappa_true", "baseline_norm_min", "baseline_norm_mean", "baseline_norm_max",
  "baseline_norm_sd", "baseline_norm_max_rel_dev", "option_user_bw_mu",
  "option_user_bw_cov", "fpca_data_type", "fpca_method_bw_mu",
  "fpca_method_bw_cov", "fpca_method_mu_cov_est", "fpca_method_xi",
  "data_warning", "truth_warning", "fit_warning", "fdapace_version",
  "dependency_md5", "r_version", "host", "slurm_job_id", "slurm_array_task_id"
)
missing_columns <- setdiff(required_columns, names(raw))
add_check("required_columns", !length(missing_columns),
          if (length(missing_columns)) paste("Missing:", paste(missing_columns, collapse = ", ")) else "All required fields are present.")
if (length(missing_columns)) {
  write_validation()
  stop("Required columns are absent; see validation_report.txt.", call. = FALSE)
}

numeric_columns <- c(
  "setting", "rep", "n", "N", "m", "seed", "config_order", "bw_multiplier",
  "kernel_bw_scale", "base_bw_mu", "base_bw_cov", "requested_bw_mu",
  "requested_bw_cov", "actual_bw_mu", "actual_bw_cov", "imse", "imse_sq",
  "eps_u", "u_cap", "Z_cap", "decayR", "beta_a_lower", "beta_a_upper",
  "beta_b_lower", "beta_b_upper", "t_grid_len", "xgrid_len", "generating_K",
  "seed_base", "pooled_time_range", "kappa_true", "baseline_norm_min",
  "baseline_norm_mean", "baseline_norm_max", "baseline_norm_sd",
  "baseline_norm_max_rel_dev", "option_user_bw_mu", "option_user_bw_cov"
)
for (name in numeric_columns) raw[[name]] <- suppressWarnings(as.numeric(raw[[name]]))
for (name in c("ok", "is_default", "user_bw_supplied", "latent")) {
  raw[[name]] <- normalize_logical(raw[[name]], name)
}
for (name in setdiff(required_columns, c(numeric_columns, "ok", "is_default", "user_bw_supplied", "latent"))) {
  raw[[name]] <- as.character(raw[[name]])
}

run_ids <- unique(raw$run_id)
expected_rows <- expected_reps * length(n_vec) * length(N_vec) * length(m_vec) *
  length(kernels) * length(bw_multipliers)
add_check("single_run_id", length(run_ids) == 1L, paste(run_ids, collapse = ","))
add_check("run_id_mode", length(run_ids) == 1L && grepl(paste0("^", mode, "_[0-9]+$"), run_ids),
          if (length(run_ids)) run_ids[1L] else "none")
add_check("setting_1_only", length(unique(raw$setting)) == 1L && unique(raw$setting) == 1,
          paste(unique(raw$setting), collapse = ","))
add_check("planned_row_count", nrow(raw) == expected_rows,
          sprintf("Observed %d; expected %d.", nrow(raw), expected_rows))

expected_grid <- expand.grid(
  run_id = if (length(run_ids)) run_ids[1L] else "", setting = 1,
  rep = seq_len(expected_reps), n = n_vec, N = N_vec, m = m_vec,
  kernel = kernels, bw_multiplier = bw_multipliers, stringsAsFactors = FALSE
)
observed_keys <- result_key(raw)
expected_keys <- paste(
  paste(expected_grid$run_id, expected_grid$setting, expected_grid$rep,
        expected_grid$n, expected_grid$N, expected_grid$m, sep = "|"),
  expected_grid$kernel, numeric_key(expected_grid$bw_multiplier), sep = "|"
)
duplicate_keys <- observed_keys[duplicated(observed_keys)]
missing_keys <- setdiff(expected_keys, observed_keys)
extra_keys <- setdiff(observed_keys, expected_keys)
add_check("unique_result_keys", !length(duplicate_keys), sprintf("Duplicates: %d.", length(unique(duplicate_keys))))
add_check("all_planned_combinations", !length(missing_keys), sprintf("Missing: %d.", length(missing_keys)))
add_check("no_unplanned_combinations", !length(extra_keys), sprintf("Unexpected: %d.", length(extra_keys)))

expected_config_id <- sprintf("%s_x%.2f", raw$kernel, raw$bw_multiplier)
add_check("configuration_ids", all(raw$config_id == expected_config_id),
          sprintf("Mismatches: %d.", sum(raw$config_id != expected_config_id)))
expected_config_order <- (kernel_order(raw$kernel) - 1L) * 5L +
  match(numeric_key(raw$bw_multiplier), numeric_key(bw_multipliers))
add_check("configuration_order", all(raw$config_order == expected_config_order),
          sprintf("Mismatches: %d.", sum(raw$config_order != expected_config_order)))
expected_default <- raw$kernel == "gauss" & close_vec(raw$bw_multiplier, 1)
add_check("default_flag", all(raw$is_default == expected_default),
          sprintf("Mismatches: %d.", sum(raw$is_default != expected_default)))

expected_seed <- 2000 + raw$rep * 1000000 + raw$n * 1000 + raw$N * 10 + raw$m
add_check("seed_formula", all(raw$seed == expected_seed), sprintf("Mismatches: %d.", sum(raw$seed != expected_seed)))
data_keys <- dataset_key(raw)
counts <- table(data_keys)
expected_datasets <- expected_reps * length(n_vec) * length(N_vec) * length(m_vec)
add_check("configurations_per_dataset", length(counts) == expected_datasets && all(counts == 15L),
          sprintf("Datasets %d/%d; configuration-count range %s.", length(counts), expected_datasets,
                  if (length(counts)) paste(range(counts), collapse = "-") else "none"))
seed_by_data <- split(raw$seed, data_keys)
one_seed <- length(seed_by_data) == expected_datasets && all(vapply(seed_by_data, function(x) length(unique(x)) == 1L, logical(1L)))
add_check("common_random_numbers", one_seed, "One seed must be shared by all 15 configurations of each dataset.")
unique_data <- raw[!duplicated(data_keys), c("rep", "n", "N", "m", "seed"), drop = FALSE]
add_check("no_seed_collisions", nrow(unique_data) == expected_datasets && !anyDuplicated(unique_data$seed),
          sprintf("Unique generated-data seeds: %d.", length(unique(unique_data$seed))))

scale_map <- c(gauss = 1, epan = sqrt(5), quar = sqrt(7))
expected_scales <- unname(scale_map[raw$kernel])
add_check("kernel_bandwidth_scales", all(close_vec(raw$kernel_bw_scale, expected_scales, 1e-12)),
          sprintf("Mismatches: %d.", sum(!close_vec(raw$kernel_bw_scale, expected_scales, 1e-12))))

base_complete <- is.finite(raw$base_bw_mu) & raw$base_bw_mu > 0 &
  is.finite(raw$base_bw_cov) & raw$base_bw_cov > 0
base_absent <- is.na(raw$base_bw_mu) & is.na(raw$base_bw_cov)
base_availability_ok <- all(base_complete | base_absent) &&
  all(vapply(split(seq_len(nrow(raw)), data_keys), function(index) {
    all(base_complete[index]) || all(base_absent[index])
  }, logical(1L)))
add_check("baseline_bandwidth_availability", base_availability_ok,
          "Each generated dataset must have one positive shared baseline pair or a wholly unavailable pair after an earlier recorded failure.")
constant_base <- base_availability_ok && all(vapply(split(seq_len(nrow(raw)), data_keys), function(index) {
  if (all(base_absent[index])) return(TRUE)
  max(raw$base_bw_mu[index]) - min(raw$base_bw_mu[index]) <= 1e-12 &&
    max(raw$base_bw_cov[index]) - min(raw$base_bw_cov[index]) <= 1e-12
}, logical(1L)))
add_check("constant_baseline_bandwidths", constant_base, "Available implicit Gaussian bandwidths must be shared within each generated dataset.")
default_fraction_ok <- all(!base_complete |
  (close_vec(raw$base_bw_mu, 0.05 * raw$pooled_time_range) &
   close_vec(raw$base_bw_cov, 0.10 * raw$pooled_time_range)))
add_check("implicit_default_bandwidth_formula", default_fraction_ok,
          "Available Gaussian baseline bandwidths must equal 0.05 and 0.10 times the pooled time range.")
expected_mu <- raw$base_bw_mu * raw$kernel_bw_scale * raw$bw_multiplier
expected_cov <- raw$base_bw_cov * raw$kernel_bw_scale * raw$bw_multiplier
request_complete <- is.finite(raw$requested_bw_mu) & raw$requested_bw_mu > 0 &
  is.finite(raw$requested_bw_cov) & raw$requested_bw_cov > 0
request_absent <- is.na(raw$requested_bw_mu) & is.na(raw$requested_bw_cov)
bandwidth_formula_ok <- all((base_complete & request_complete &
  close_vec(raw$requested_bw_mu, expected_mu) & close_vec(raw$requested_bw_cov, expected_cov)) |
  (base_absent & request_absent))
add_check("requested_bandwidth_formula", bandwidth_formula_ok,
          "Available requests must equal base h times kernel scale times multiplier; unavailable references must have no request.")

fixed_ok <-
  all(raw$mode == mode) && all(raw$setting_label == "Setting 1") &&
  all(!raw$latent) && all(raw$methodXi == "CE") && all(raw$methodMuCovEst == "smooth") &&
  all(close_vec(raw$eps_u, 0.05)) && all(close_vec(raw$u_cap, 0.95)) &&
  all(close_vec(raw$Z_cap, 8)) && all(close_vec(raw$decayR, 1)) &&
  all(close_vec(raw$beta_a_lower, 0.6)) && all(close_vec(raw$beta_a_upper, 3)) &&
  all(close_vec(raw$beta_b_lower, 0.6)) && all(close_vec(raw$beta_b_upper, 3)) &&
  all(raw$t_design == "random") && all(raw$t_grid_len == 51) &&
  all(raw$xgrid_len == 128) && all(raw$generating_K == 50) && all(raw$seed_base == 2000) &&
  all(raw$fdapace_version == "0.6.0") &&
  all(raw$dependency_md5 == "4beac3d60a2b3e8850037ea7993abc05") &&
  length(unique(raw$r_version)) == 1L && all(raw$r_version == "4.5.0") &&
  length(unique(raw$host)) == 1L && all(nzchar(raw$host))
add_check("row_scientific_provenance", fixed_ok, "All fixed DGP, estimator, package, and dependency fields must match the fixed design.")

slurm_ids <- unique(raw$slurm_job_id)
slurm_ok <- length(slurm_ids) == 1L && nzchar(slurm_ids) &&
  length(run_ids) == 1L && identical(run_ids, paste0(mode, "_", slurm_ids)) &&
  all(is.na(raw$slurm_array_task_id) | !nzchar(raw$slurm_array_task_id))
add_check("slurm_identity", slurm_ok, "Row Slurm identity must be constant, nonempty, non-array, and agree with run_id.")

config_ok <- FALSE
if (config_readable) {
  cfg <- config_bundle$config
  grid <- config_bundle$tuning_grid
  expected_grid_kernel <- rep(kernels, each = length(bw_multipliers))
  expected_grid_multiplier <- rep(bw_multipliers, times = length(kernels))
  expected_grid_id <- sprintf("%s_x%.2f", expected_grid_kernel, expected_grid_multiplier)
  config_ok <-
    identical(as.character(cfg$run_id), run_ids) && identical(cfg$mode, mode) &&
    identical(as.integer(cfg$setting_id), 1L) &&
    identical(as.integer(cfg$reps), as.integer(expected_reps)) &&
    identical(as.integer(cfg$n_vec), as.integer(n_vec)) &&
    identical(as.integer(cfg$N_vec), as.integer(N_vec)) &&
    identical(as.integer(cfg$m_vec), as.integer(m_vec)) &&
    identical(as.character(cfg$kernels), kernels) &&
    same_numeric(cfg$bw_multipliers, bw_multipliers) &&
    identical(as.integer(cfg$seed_base), 2000L) && identical(cfg$latent, FALSE) &&
    identical(as.character(cfg$methodXi), "CE") &&
    identical(as.character(cfg$methodMuCovEst), "smooth") && same_numeric(cfg$decayR, 1) &&
    same_numeric(cfg$beta_a, c(0.6, 3)) && same_numeric(cfg$beta_b, c(0.6, 3)) &&
    identical(as.character(cfg$t_design), "random") && same_numeric(cfg$eps_u, 0.05) &&
    same_numeric(cfg$u_cap, 0.95) && same_numeric(cfg$Z_cap, 8) &&
    identical(as.integer(cfg$t_grid_len), 51L) && identical(as.integer(cfg$xgrid_len), 128L) &&
    identical(as.integer(cfg$generating_K), 50L) &&
    identical(as.character(cfg$fdapace_version), "0.6.0") &&
    identical(as.character(cfg$dependency_md5), "4beac3d60a2b3e8850037ea7993abc05") &&
    identical(as.character(cfg$slurm_job_id), slurm_ids) &&
    identical(as.integer(cfg$planned_rows), as.integer(expected_rows)) &&
    isTRUE(cfg$explicit_default_option_equivalence) && nrow(grid) == 15L &&
    identical(as.character(grid$kernel), expected_grid_kernel) &&
    same_numeric(grid$bw_multiplier, expected_grid_multiplier) &&
    identical(as.character(grid$config_id), expected_grid_id) &&
    identical(as.integer(grid$config_order), 1:15) &&
    identical(as.logical(grid$is_default), expected_grid_kernel == "gauss" & expected_grid_multiplier == 1) &&
    all(close_vec(grid$kernel_bw_scale, unname(scale_map[grid$kernel]), 1e-12))
}
add_check("run_config_scientific_provenance", config_ok, "run_config.rds must encode the exact fixed design and kernel-scale map.")

metadata_checksum_ok <- if (combined_ok) {
  identical(as.character(metadata$raw_results_rds_md5), unname(tools::md5sum(combined_path)))
} else !combined_present && rep_read_ok
metadata_ok <- metadata_readable && identical(as.integer(metadata$planned_rows), as.integer(expected_rows)) &&
  identical(as.integer(metadata$observed_rows), as.integer(nrow(raw))) &&
  identical(as.integer(metadata$rep_files), as.integer(expected_reps)) &&
  identical(as.integer(metadata$unique_seeds), as.integer(expected_datasets)) &&
  isTRUE(metadata$structural_validation_passed) && metadata_checksum_ok
add_check("run_metadata_consistency", metadata_ok, "Driver metadata, row counts, seeds, replication files, and raw checksum must agree.")

default_rows <- raw$is_default
explicit_rows <- !default_rows
available_default <- default_rows & base_complete
unavailable_default <- default_rows & base_absent
available_explicit <- explicit_rows & base_complete
unavailable_explicit <- explicit_rows & base_absent
default_contract <- all(raw$bandwidth_mode[available_default] == "implicit_default") &&
  all(raw$bandwidth_mode[unavailable_default] == "unavailable") &&
  all(!raw$user_bw_supplied[default_rows]) &&
  all(raw$actual_kernel[default_rows & raw$ok] == "gauss") &&
  all(raw$fpca_method_bw_mu[default_rows & raw$ok] == "Default") &&
  all(raw$fpca_method_bw_cov[default_rows & raw$ok] == "Default")
explicit_contract <- all(raw$bandwidth_mode[available_explicit] == "explicit_second_moment_scaled") &&
  all(raw$bandwidth_mode[unavailable_explicit] == "unavailable") &&
  all(raw$user_bw_supplied[available_explicit]) &&
  all(!raw$user_bw_supplied[unavailable_explicit]) &&
  all(raw$fpca_method_bw_mu[explicit_rows & raw$ok] == "Default") &&
  all(raw$fpca_method_bw_cov[explicit_rows & raw$ok] == "Default")
add_check("default_execution_contract", default_contract, "The reference must be the untouched implicit Gaussian fit.")
add_check("explicit_execution_contract", explicit_contract, "All nonreference rows must use explicit scaled bandwidths.")

successful <- raw$ok
actual_ok <- all(!successful | (
  raw$actual_kernel == raw$kernel & close_vec(raw$actual_bw_mu, raw$requested_bw_mu, 1e-7) &
    close_vec(raw$actual_bw_cov, raw$requested_bw_cov, 1e-7) &
    close_vec(raw$option_user_bw_mu, raw$requested_bw_mu, 1e-7) &
    close_vec(raw$option_user_bw_cov, raw$requested_bw_cov, 1e-7) &
    raw$fpca_data_type == "Sparse" & raw$fpca_method_mu_cov_est == "smooth" &
    raw$fpca_method_xi == "CE"
))
add_check("returned_fit_contract", actual_ok, "Successful fits must return the requested kernel/bandwidth and fixed FPCA options.")

baseline_columns <- c("baseline_norm_min", "baseline_norm_mean", "baseline_norm_max",
                      "baseline_norm_sd", "baseline_norm_max_rel_dev")
baseline_finite <- apply(raw[baseline_columns], 1L, function(x) all(is.finite(as.numeric(x))))
simulation_available <- is.finite(raw$kappa_true)
baseline_absent <- apply(raw[baseline_columns], 1L, function(x) all(is.na(as.numeric(x))))
baseline_order_ok <- simulation_available & baseline_finite & raw$kappa_true > 0 & raw$baseline_norm_min >= 0 &
  raw$baseline_norm_sd >= 0 & raw$baseline_norm_max_rel_dev >= 0 &
  raw$baseline_norm_min <= raw$baseline_norm_mean & raw$baseline_norm_mean <= raw$baseline_norm_max
expected_rel_dev <- pmax(abs(raw$baseline_norm_min - raw$kappa_true),
                         abs(raw$baseline_norm_max - raw$kappa_true)) / abs(raw$kappa_true)
baseline_formula_ok <- close_vec(raw$baseline_norm_max_rel_dev, expected_rel_dev)
baseline_by_data <- unique(raw[c("run_id", "setting", "rep", "n", "N", "m", "seed",
                                 "kappa_true", baseline_columns)])
baseline_ok <- all((simulation_available & baseline_order_ok & baseline_formula_ok) |
                   (!simulation_available & baseline_absent)) &&
  nrow(baseline_by_data) == expected_datasets
add_check("baseline_norm_diagnostics", baseline_ok,
          sprintf("Unique dataset diagnostics: %d; expected %d.", nrow(baseline_by_data), expected_datasets))

time_available <- is.finite(raw$pooled_time_range) & raw$pooled_time_range > 0
stage_allowed <- raw$stage %in% c("simulation", "truth", "fit_default", "reference_unavailable",
                                  "fit", "evaluation", "complete")
stage_ok <- all(stage_allowed) && all(raw$stage[raw$ok] == "complete") &&
  all(time_available == (raw$stage != "simulation")) &&
  all(simulation_available == (raw$stage != "simulation"))
add_check("stage_and_time_provenance", stage_ok,
          "Stages must follow the driver contract, with time/baseline diagnostics absent only for simulation failures.")

status_ok <- all(raw$status == ifelse(raw$ok, "success", "failure")) &&
  all(nzchar(raw$err[!raw$ok])) && all(!nzchar(raw$err[raw$ok]))
add_check("row_status_and_errors", status_ok, "Success/failure labels and error messages must agree.")
finite_metrics <- is.finite(raw$imse[successful]) & is.finite(raw$imse_sq[successful])
nonnegative_metrics <- raw$imse[successful] >= 0 & raw$imse_sq[successful] >= 0
add_check("finite_successful_metrics", all(finite_metrics), sprintf("Invalid rows: %d.", sum(!finite_metrics)))
add_check("nonnegative_successful_metrics", all(nonnegative_metrics), sprintf("Negative rows: %d.", sum(!nonnegative_metrics)))
add_check("failed_metrics_missing", all(is.na(raw$imse[!successful]) & is.na(raw$imse_sq[!successful])),
          "Failed fits must not carry fabricated performance values.")
failure_count <- sum(!successful)
expected_compact_window_messages <- paste0(
  "Explicit sensitivity fit failed: ",
  c("No enough points in local window, please increase bandwidth.",
    "No enough points in local window, please increase bandwidth using userBwCov.")
)
expected_compact_window_failure <- !successful & raw$stage == "fit" &
  raw$kernel %in% c("epan", "quar") &
  raw$err %in% expected_compact_window_messages
unexpected_failure_count <- sum(!successful & !expected_compact_window_failure)
add_check("no_unexpected_pipeline_failures", unexpected_failure_count == 0L,
          sprintf("Unexpected simulation/truth/default/evaluation/fit failures: %d.", unexpected_failure_count))
add_check("fit_failures_accounted", TRUE,
          sprintf("Recorded compact-kernel local-window failures: %d; treated as a numerical outcome, not a structural omission.",
                  sum(expected_compact_window_failure)))
numerical_status <- if (failure_count == 0L) "COMPLETE" else "COMPLETE_WITH_FIT_FAILURES"

## Structural/provenance checks must pass before numerical summaries are created.
structural_pass <- all(vapply(checks, function(x) isTRUE(x$passed), logical(1L)))

failure_df <- raw[!successful, c("run_id", "rep", "n", "N", "m", "seed", "kernel",
                                 "bw_multiplier", "stage", "err"), drop = FALSE]
failure_summary <- if (nrow(failure_df)) {
  tmp <- aggregate(rep(1L, nrow(failure_df)), list(error = failure_df$err), sum)
  names(tmp)[2L] <- "count"
  tmp[order(-tmp$count, tmp$error), , drop = FALSE]
} else data.frame(error = character(), count = integer(), stringsAsFactors = FALSE)
atomic_write_csv(failure_df, file.path(outdir, "failed_rows.csv"))
atomic_write_csv(failure_summary, file.path(outdir, "failure_summary.csv"))

if (!structural_pass) {
  overall <- write_validation()
  stop("Structural/provenance validation failed; numerical summaries were not generated. See ",
       file.path(outdir, "validation_report.txt"), call. = FALSE)
}

atomic_write_lines(c(
  paste0("status=", numerical_status),
  paste0("planned_rows=", nrow(raw)),
  paste0("successful_rows=", sum(successful)),
  paste0("failed_rows=", failure_count)
), file.path(outdir, "numerical_status.txt"))

raw <- sort_results(raw, include_rep = TRUE)
atomic_write_rds(raw, file.path(outdir, "checked_raw_results.rds"))
atomic_write_csv(raw, file.path(outdir, "checked_raw_results.csv"))

baseline_by_data <- baseline_by_data[order(baseline_by_data$rep, baseline_by_data$n,
                                           baseline_by_data$N, baseline_by_data$m), , drop = FALSE]
atomic_write_csv(baseline_by_data, file.path(outdir, "baseline_norm_by_dataset.csv"))
baseline_summary <- group_apply(baseline_by_data, c("run_id", "setting", "n", "N", "m"), function(d) {
  data.frame(
    datasets = nrow(d),
    mean_baseline_norm = mean(d$baseline_norm_mean),
    mcse_dataset_mean_baseline_norm = mcse(d$baseline_norm_mean),
    mean_max_relative_deviation = mean(d$baseline_norm_max_rel_dev),
    max_relative_deviation = max(d$baseline_norm_max_rel_dev),
    stringsAsFactors = FALSE
  )
})
atomic_write_csv(baseline_summary, file.path(outdir, "baseline_norm_summary.csv"))

warning_columns <- c("data_warning", "truth_warning", "fit_warning")
warning_any <- apply(raw[warning_columns], 1L, function(x) any(nzchar(as.character(x))))
warning_rows <- raw[warning_any, c("run_id", "rep", "n", "N", "m", "kernel",
                                  "bw_multiplier", warning_columns), drop = FALSE]
atomic_write_csv(warning_rows, file.path(outdir, "warning_rows.csv"))

summary_columns <- c("run_id", "setting", "kernel", "bw_multiplier", "config_id", "n", "N", "m")
summary_by_cell <- group_apply(raw, summary_columns, function(d) {
  ok <- d$ok & is.finite(d$imse) & is.finite(d$imse_sq)
  data.frame(
    reps_planned = nrow(d), reps_ok = sum(ok), failures = sum(!ok),
    mean_imse = safe_mean(d$imse[ok]), sd_imse = safe_sd(d$imse[ok]), mcse_imse = mcse(d$imse[ok]),
    mean_imse_sq = safe_mean(d$imse_sq[ok]), sd_imse_sq = safe_sd(d$imse_sq[ok]),
    mcse_imse_sq = mcse(d$imse_sq[ok]),
    mean_base_bw_mu = safe_mean(d$base_bw_mu), mean_base_bw_cov = safe_mean(d$base_bw_cov),
    mean_requested_bw_mu = safe_mean(d$requested_bw_mu), mean_requested_bw_cov = safe_mean(d$requested_bw_cov),
    mean_actual_bw_mu = safe_mean(d$actual_bw_mu[ok]), mean_actual_bw_cov = safe_mean(d$actual_bw_cov[ok]),
    stringsAsFactors = FALSE
  )
})
summary_by_cell <- sort_results(summary_by_cell)
atomic_write_csv(summary_by_cell, file.path(outdir, "summary_by_cell.csv"))

reference <- raw[raw$is_default, c("run_id", "setting", "rep", "n", "N", "m", "seed",
                                    "ok", "imse", "imse_sq"), drop = FALSE]
names(reference)[names(reference) == "ok"] <- "default_ok"
names(reference)[names(reference) == "imse"] <- "default_imse"
names(reference)[names(reference) == "imse_sq"] <- "default_imse_sq"
pair_by <- c("run_id", "setting", "rep", "n", "N", "m", "seed")
paired <- merge(raw, reference, by = pair_by, all.x = TRUE, sort = FALSE)
paired$pair_ok <- paired$ok & paired$default_ok & is.finite(paired$imse) &
  is.finite(paired$default_imse) & paired$default_imse > 0 & is.finite(paired$imse_sq) &
  is.finite(paired$default_imse_sq) & paired$default_imse_sq > 0
paired$paired_diff_imse_sq <- ifelse(paired$pair_ok, paired$imse_sq - paired$default_imse_sq, NA_real_)
paired$paired_pct_imse_sq <- ifelse(paired$pair_ok, 100 * paired$paired_diff_imse_sq / paired$default_imse_sq, NA_real_)
paired$ratio_imse_sq <- ifelse(paired$pair_ok, paired$imse_sq / paired$default_imse_sq, NA_real_)
paired$paired_diff_imse <- ifelse(paired$pair_ok, paired$imse - paired$default_imse, NA_real_)
paired$paired_pct_imse <- ifelse(paired$pair_ok, 100 * paired$paired_diff_imse / paired$default_imse, NA_real_)
paired$ratio_imse <- ifelse(paired$pair_ok, paired$imse / paired$default_imse, NA_real_)
paired <- sort_results(paired, include_rep = TRUE)
atomic_write_rds(paired, file.path(outdir, "paired_replication_results.rds"))
atomic_write_csv(paired, file.path(outdir, "paired_replication_results.csv"))

paired_summary <- group_apply(paired, summary_columns, function(d) {
  ok <- d$pair_ok
  config_mean_sq <- safe_mean(d$imse_sq[ok])
  default_mean_sq <- safe_mean(d$default_imse_sq[ok])
  config_mean <- safe_mean(d$imse[ok])
  default_mean <- safe_mean(d$default_imse[ok])
  data.frame(
    pairs_planned = nrow(d), pairs_ok = sum(ok), failures = nrow(d) - sum(ok),
    mean_imse_sq = config_mean_sq, mcse_imse_sq = mcse(d$imse_sq[ok]),
    default_mean_imse_sq = default_mean_sq,
    relative_change_pct = if (is.finite(default_mean_sq) && default_mean_sq > 0)
      100 * (config_mean_sq / default_mean_sq - 1) else NA_real_,
    paired_mean_diff_imse_sq = safe_mean(d$paired_diff_imse_sq[ok]),
    paired_mcse_diff_imse_sq = mcse(d$paired_diff_imse_sq[ok]),
    paired_mean_pct_imse_sq = safe_mean(d$paired_pct_imse_sq[ok]),
    paired_mcse_pct_imse_sq = mcse(d$paired_pct_imse_sq[ok]),
    mean_ratio_imse_sq = safe_mean(d$ratio_imse_sq[ok]),
    mean_imse = config_mean, mcse_imse = mcse(d$imse[ok]), default_mean_imse = default_mean,
    relative_change_imse_pct = if (is.finite(default_mean) && default_mean > 0)
      100 * (config_mean / default_mean - 1) else NA_real_,
    paired_mean_diff_imse = safe_mean(d$paired_diff_imse[ok]),
    paired_mcse_diff_imse = mcse(d$paired_diff_imse[ok]),
    paired_mean_pct_imse = safe_mean(d$paired_pct_imse[ok]),
    paired_mcse_pct_imse = mcse(d$paired_pct_imse[ok]),
    stringsAsFactors = FALSE
  )
})
paired_summary <- sort_results(paired_summary)
atomic_write_csv(paired_summary, file.path(outdir, "paired_summary_by_cell.csv"))

across_columns <- c("run_id", "setting", "kernel", "bw_multiplier", "config_id")
summary_across_cells <- group_apply(paired_summary, across_columns, function(d) {
  rel <- d$relative_change_pct
  valid <- is.finite(rel)
  abs_rel <- abs(rel[valid])
  worst <- if (any(valid)) which(valid)[which.max(abs_rel)] else NA_integer_
  data.frame(
    cells_planned = nrow(d), cells_ok = sum(valid),
    mean_relative_change_pct = safe_mean(rel),
    median_relative_change_pct = if (any(valid)) stats::median(rel[valid]) else NA_real_,
    mean_abs_relative_change_pct = safe_mean(abs(rel)),
    median_abs_relative_change_pct = if (any(valid)) stats::median(abs_rel) else NA_real_,
    max_abs_relative_change_pct = if (any(valid)) max(abs_rel) else NA_real_,
    min_relative_change_pct = if (any(valid)) min(rel[valid]) else NA_real_,
    max_relative_change_pct = if (any(valid)) max(rel[valid]) else NA_real_,
    cells_small_le_5pct = sum(valid & abs(rel) <= 5),
    cells_moderate_gt5_le10pct = sum(valid & abs(rel) > 5 & abs(rel) <= 10),
    cells_material_gt10pct = sum(valid & abs(rel) > 10),
    worst_n = if (is.na(worst)) NA_integer_ else d$n[worst],
    worst_N = if (is.na(worst)) NA_integer_ else d$N[worst],
    worst_m = if (is.na(worst)) NA_integer_ else d$m[worst],
    worst_relative_change_pct = if (is.na(worst)) NA_real_ else rel[worst],
    total_failed_pairs = sum(d$pairs_planned - d$pairs_ok),
    stringsAsFactors = FALSE
  )
})
summary_across_cells <- sort_results(summary_across_cells)
atomic_write_csv(summary_across_cells, file.path(outdir, "summary_across_cells.csv"))

representative_abs <- summary_by_cell[summary_by_cell$n == 100 & summary_by_cell$N == 5 &
                                        summary_by_cell$m == 200, , drop = FALSE]
representative_pair <- paired_summary[paired_summary$n == 100 & paired_summary$N == 5 &
                                        paired_summary$m == 200, , drop = FALSE]
if (mode == "full") {
  representative_ok <- nrow(representative_abs) == 15L && nrow(representative_pair) == 15L
  add_check("representative_table_cell", representative_ok,
            sprintf("Observed %d absolute and %d paired configurations; expected 15 each.",
                    nrow(representative_abs), nrow(representative_pair)))
}

failures_by_config <- group_apply(raw, c("run_id", "setting", "kernel", "bw_multiplier", "config_id"),
                                  function(d) data.frame(
                                    fits_planned = nrow(d),
                                    fits_ok = sum(d$ok),
                                    failures = sum(!d$ok),
                                    failure_pct = 100 * sum(!d$ok) / nrow(d)
                                  ))
failures_by_config <- sort_results(failures_by_config)
atomic_write_csv(failures_by_config, file.path(outdir, "failure_summary_by_config.csv"))
if (write_table) {
  if (mode != "full") {
    add_check("latex_table_generated", FALSE, "A numerical table is generated only for the full design.")
  } else {
    table_script <- file.path(get_script_dir(), "generate_table.R")
    if (!file.exists(table_script)) stop("Cannot find table renderer: ", table_script)
    table_env <- new.env(parent = globalenv())
    sys.source(table_script, envir = table_env)
    table_env$write_sensitivity_table(
      representative_abs, representative_pair, summary_across_cells,
      failures_by_config, table_out
    )
    add_check("latex_table_generated", file.exists(table_out), normalizePath(table_out, mustWork = TRUE))
  }
}

overall_pass <- write_validation()
cat("Structural validation: ", if (overall_pass) "PASS" else "FAIL", "\n", sep = "")
cat("Numerical completion: ", numerical_status, "\n", sep = "")
cat("Rows: ", nrow(raw), "\n", sep = "")
cat("Failed fits: ", failure_count, "\n", sep = "")
cat("Aggregation outputs: ", outdir, "\n", sep = "")
if (write_table && file.exists(table_out)) cat("LaTeX table: ", table_out, "\n", sep = "")

if (strict && !overall_pass) {
  stop("Sensitivity validation failed; see ", file.path(outdir, "validation_report.txt"), call. = FALSE)
}
