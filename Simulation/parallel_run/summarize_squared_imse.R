#!/usr/bin/env Rscript

# Figure 3 and Table 1 use mean_i integral_t W1_i(t)^2 dt.
# Approved publication comparison: proposed latent=TRUE/methodXi=CE
# (estimator (14)) versus LQD and Tangent with methodXi=IN.
# Source this file to use read_squared_results() without running summaries.
# Rscript summarize_squared_imse.R writes only to the publication output directory.

read_squared_results <- function(prop_dir, base_dir, expected_reps = 200L, setting = 1L) {
  if (length(expected_reps) != 1L || !is.numeric(expected_reps) ||
      !is.finite(expected_reps) || expected_reps < 1 ||
      expected_reps != floor(expected_reps)) {
    stop("expected_reps must be a positive integer")
  }
  if (length(setting) != 1L || is.na(setting) || !setting %in% c(1L, 2L)) {
    stop("setting must be 1 or 2")
  }
  common_metadata <- list(
    reps = expected_reps, seed_base = 2000,
    t_design = "random", decayR = 1,
    n_vec = c(20, 50, 100, 200, 400), N_vec = c(3, 5, 10, 20), m_vec = 200,
    beta_a = if (setting == 1L) c(.6, 3) else c(3, 4),
    beta_b = if (setting == 1L) c(.6, 3) else c(1, 2),
    t_grid_out = seq(0, 1, length.out = 51)
  )
  check_metadata <- function(directory, filename, specific, label) {
    metadata_path <- file.path(directory, filename)
    if (!file.exists(metadata_path)) stop("Missing ", label, " metadata: ", metadata_path)
    metadata <- new.env(parent = emptyenv())
    load(metadata_path, envir = metadata)
    expected_metadata <- c(common_metadata, specific)
    for (field in names(expected_metadata)) {
      if (!exists(field, metadata, inherits = FALSE)) {
        stop(label, " missing metadata field ", field, ", setting ", setting)
      }
      expected_value <- expected_metadata[[field]]
      actual <- if (exists(field, metadata, inherits = FALSE)) get(field, metadata) else NULL
      matches <- if (is.numeric(expected_value)) {
        is.numeric(actual) && isTRUE(all.equal(as.numeric(actual),
                                              as.numeric(expected_value), tolerance = 1e-12))
      } else identical(actual, expected_value)
      if (!matches) stop(label, " metadata mismatch for ", field, ", setting ", setting)
    }
    if (!exists("raw_df", metadata, inherits = FALSE) ||
        !identical(get("raw_df", metadata), readRDS(file.path(directory, "raw_results.rds")))) {
      stop(label, " raw_results.rds does not match its saved metadata bundle")
    }
  }
  check_metadata(prop_dir, "results_prop.RData",
                 list(latent = TRUE, methodXi = "CE", eps_u = .05, u_cap = .95, Z_cap = 8),
                 "Proposed latent CE")
  check_metadata(base_dir, "results_baselines.RData",
                 list(methodXi = "IN", Ku = 10, basis_type = "bspline", bs_degree = 3,
                      smooth_deriv = FALSE, spar = NULL), "IN benchmarks")
  keys <- c("rep", "n", "N", "m")
  expected <- expand.grid(rep = seq_len(expected_reps),
                          n = c(20, 50, 100, 200, 400),
                          N = c(3, 5, 10, 20), m = 200)
  make_key <- function(d) do.call(paste, c(d[keys], sep = "/"))
  expected_keys <- make_key(expected)

  read_checked <- function(directory, metrics, label) {
    path <- file.path(directory, "raw_results.rds")
    if (!file.exists(path)) stop(label, ": missing ", path)
    raw <- readRDS(path)
    if (!is.data.frame(raw)) stop(label, ": raw_results.rds must be a data frame")
    if (anyDuplicated(names(raw))) stop(label, ": duplicate column names")
    required <- c(keys, "ok", metrics)
    missing <- setdiff(required, names(raw))
    if (length(missing)) stop(label, ": missing columns: ", paste(missing, collapse = ", "))
    for (key in keys) {
      values <- raw[[key]]
      if (!is.numeric(values) || any(!is.finite(values)) ||
          any(values != floor(values))) {
        stop(label, ": ", key, " must contain finite integers")
      }
    }
    actual_keys <- make_key(raw)
    if (anyDuplicated(actual_keys)) {
      stop(label, ": duplicate (rep,n,N,m) key: ", actual_keys[anyDuplicated(actual_keys)])
    }
    missing_keys <- setdiff(expected_keys, actual_keys)
    extra_keys <- setdiff(actual_keys, expected_keys)
    if (length(missing_keys) || length(extra_keys)) {
      stop(label, ": incomplete or unexpected grid; missing ", length(missing_keys),
           " keys, extra ", length(extra_keys), " keys",
           if (length(missing_keys)) paste0("; first missing: ", missing_keys[1]) else "",
           if (length(extra_keys)) paste0("; first extra: ", extra_keys[1]) else "")
    }
    if (!is.logical(raw$ok) || anyNA(raw$ok) || !all(raw$ok)) {
      stop(label, ": all planned rows must have ok=TRUE; failures are not dropped")
    }
    for (metric in metrics) {
      values <- raw[[metric]]
      if (!is.numeric(values) || any(!is.finite(values)) || any(values < 0)) {
        stop(label, ": ", metric, " must be finite and nonnegative in every planned row")
      }
    }
    raw[match(expected_keys, actual_keys), c(keys, metrics), drop = FALSE]
  }

  proposed <- read_checked(prop_dir, "imse_sq", "Proposed")
  baseline <- read_checked(base_dir, c("imse_sq_lqd", "imse_sq_tangent"), "Baselines")
  # Both tables have already been matched to the identical, complete grid.
  merged <- proposed[keys]
  merged$imse_sq_prop <- proposed$imse_sq
  merged$imse_sq_lqd <- baseline$imse_sq_lqd
  merged$imse_sq_tangent <- baseline$imse_sq_tangent
  merged <- merged[order(merged$rep, merged$n, merged$N, merged$m), , drop = FALSE]
  rownames(merged) <- NULL
  merged
}

summarize_squared_imse <- function(results_dir, output_dir, expected_reps = 200L) {
  settings <- c("setting1", "setting2")
  # Validate both settings before writing any output.
  records <- lapply(seq_along(settings), function(index) {
    read_squared_results(file.path(results_dir, "proposed_latent_CE", settings[index]),
                         file.path(results_dir, "benchmarks_IN", settings[index]),
                         expected_reps, setting = index)
  })
  names(records) <- settings
  methods <- c(Proposed = "imse_sq_prop", LQD = "imse_sq_lqd",
               Tangent = "imse_sq_tangent")
  summary_rows <- list()
  improvement_rows <- list()
  for (setting_index in seq_along(settings)) {
    data <- records[[setting_index]]
    for (n0 in c(20, 50, 100, 200, 400)) {
      for (N0 in c(3, 5, 10, 20)) {
        d <- data[data$n == n0 & data$N == N0 & data$m == 200, , drop = FALSE]
        means <- vapply(methods, function(column) mean(d[[column]]), numeric(1))
        if (any(means[c("LQD", "Tangent")] <= 0)) {
          stop(settings[setting_index], ": percentage improvement is undefined for a zero benchmark mean")
        }
        for (method in names(methods)) {
          values <- d[[methods[[method]]]]
          summary_rows[[length(summary_rows) + 1L]] <- data.frame(
            setting = setting_index, n = n0, N = N0, m = 200,
            method = method, reps = nrow(d),
            mean_imse_sq = mean(values),
            mcse_imse_sq = stats::sd(values) / sqrt(nrow(d)),
            median_imse_sq = stats::median(values),
            sd_imse_sq = stats::sd(values), iqr_imse_sq = stats::IQR(values),
            stringsAsFactors = FALSE
          )
        }
        improvement_rows[[length(improvement_rows) + 1L]] <- data.frame(
          setting = setting_index, n = n0, N = N0, m = 200, reps = nrow(d),
          mean_imse_sq_prop = unname(means["Proposed"]),
          mean_imse_sq_lqd = unname(means["LQD"]),
          mean_imse_sq_tangent = unname(means["Tangent"]),
          improvement_lqd_pct = unname(100 * (1 - means["Proposed"] / means["LQD"])),
          improvement_tangent_pct = unname(100 * (1 - means["Proposed"] / means["Tangent"])),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  summary <- do.call(rbind, summary_rows)
  improvements <- do.call(rbind, improvement_rows)
  rownames(summary) <- rownames(improvements) <- NULL

  table_rows <- character()
  for (setting_index in seq_along(settings)) {
    if (setting_index > 1L) table_rows <- c(table_rows, "\\midrule")
    table_rows <- c(table_rows, sprintf("\\multirow{5}{*}{Setting %d}", setting_index))
    for (n0 in c(20, 50, 100, 200, 400)) {
      d <- improvements[improvements$setting == setting_index & improvements$n == n0, ]
      # The publication table matches Figure 3's N=3,5,10 panels.
      # N=20 remains in the complete summaries below.
      d <- d[match(c(3, 5, 10), d$N), ]
      values <- c(d$improvement_lqd_pct, d$improvement_tangent_pct)
      numbers <- sprintf("\\hang{%.2f}", values)
      table_rows <- c(table_rows, paste0("& ", n0, " & ",
                                        paste(numbers, collapse = " & "), " \\\\"))
    }
  }

  validation <- c(
    "PASS: squared-W1 IMSE inputs validated without dropping rows.",
    sprintf("Each method/setting has exactly %d unique rows: reps 1:%d x n={20,50,100,200,400} x N={3,5,10,20} x m=200.",
            expected_reps * 20L, expected_reps),
    "Both proposed and baseline tables contain every planned key and no extra keys.",
    "All ok values are TRUE; all squared-error values are finite and nonnegative.",
    "Criterion: mean_i integral_t [integral_u |estimated transport - true transport| du]^2 dt.",
    "Proposed input: imse_sq; baseline inputs: imse_sq_lqd and imse_sq_tangent.",
    "Proposed metadata verifies estimator (14): latent=TRUE, methodXi=CE, eps_u=.05, u_cap=.95, Z_cap=8.",
    "Benchmark metadata verifies methodXi=IN for both LQD and Tangent, with Ku=10 and the original spline options.",
    "Both methods' metadata verify original seeds/grid and intended Beta bounds in both settings.",
    "Each raw_results.rds is identical to raw_df in its verified metadata bundle.",
    "Publication Table 1 displays N={3,5,10}; complete summaries retain all N={3,5,10,20} results.",
    "Table percentages: 100 * (1 - Monte Carlo mean proposed IMSE / Monte Carlo mean benchmark IMSE).",
    "All benchmark means are positive.",
    "Summary statistics include Monte Carlo mean, MCSE, median, SD and IQR; SD/IQR describe replication variability.",
    "These summaries do not change any estimator, fitted trajectory, or simulation setting."
  )
  input_paths <- unlist(lapply(settings, function(setting) {
    file.path(results_dir, c("proposed_latent_CE", "benchmarks_IN"), setting, "raw_results.rds")
  }))
  input_paths <- c(input_paths,
                  file.path(results_dir, "proposed_latent_CE", settings, "results_prop.RData"),
                  file.path(results_dir, "benchmarks_IN", settings, "results_baselines.RData"))
  validation <- c(validation, "Input MD5 checksums:",
                  paste(unname(tools::md5sum(input_paths)), normalizePath(input_paths)))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(summary, file.path(output_dir, "summary_by_cell.csv"), row.names = FALSE)
  utils::write.csv(improvements, file.path(output_dir, "table1_improvement.csv"), row.names = FALSE)
  utils::write.csv(improvements[improvements$N %in% c(3, 5, 10), ],
                   file.path(output_dir, "table1_display.csv"), row.names = FALSE)
  writeLines(table_rows, file.path(output_dir, "table1_rows.tex"))
  writeLines(validation, file.path(output_dir, "validation_report.txt"))
  invisible(list(summary = summary, improvements = improvements, validation = validation))
}

if (sys.nframe() == 0L) {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) != 1L) stop("Cannot resolve this script's directory")
  script_dir <- dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg), fixed = TRUE)))
  args <- commandArgs(trailingOnly = TRUE)
  if (any(!grepl("=", args, fixed = TRUE))) stop("Arguments must use key=value.")
  kv <- setNames(as.list(sub("^[^=]*=", "", args)), sub("=.*$", "", args))
  if (length(setdiff(names(kv), c("results_dir", "outdir")))) stop("Unknown argument.")
  results_dir <- if (is.null(kv$results_dir)) file.path(script_dir, "results", "paper_comparison") else kv$results_dir
  output_dir <- if (is.null(kv$outdir)) file.path(results_dir, "summary") else kv$outdir
  summarize_squared_imse(results_dir, output_dir)
  cat("Wrote squared-W1 IMSE summaries to", normalizePath(output_dir), "\n")
}
