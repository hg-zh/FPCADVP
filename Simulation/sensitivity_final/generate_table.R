#!/usr/bin/env Rscript

## Render the supplement table from validated Setting-1 summary files.

options(stringsAsFactors = FALSE)

parse_kv_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) {
      stop("Every argument must have key=value form: ", arg)
    }
    key <- sub("=.*$", "", arg)
    value <- sub("^[^=]*=", "", arg)
    if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", key) || !is.null(out[[key]])) {
      stop("Invalid or duplicated argument key: ", key)
    }
    out[[key]] <- value
  }
  out
}

numeric_key <- function(x) sprintf("%.8f", as.numeric(x))

planned_grid <- function() {
  do.call(rbind, lapply(c("gauss", "epan", "quar"), function(kernel) {
    data.frame(
      kernel = kernel,
      bw_multiplier = c(0.50, 0.75, 1.00, 1.25, 1.50),
      stringsAsFactors = FALSE
    )
  }))
}

config_key <- function(x) paste(x$kernel, numeric_key(x$bw_multiplier), sep = "|")

require_columns <- function(x, columns, source) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(source, " is missing required column(s): ", paste(missing, collapse = ", "))
  }
}

require_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

integer_like <- function(x) {
  is.finite(x) & abs(x - round(x)) <= 1e-12
}

atomic_write_lines <- function(lines, path) {
  if (file.exists(path)) stop("Refusing to overwrite existing table: ", path)
  parent <- dirname(path)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create output directory: ", parent)
  }
  temporary <- tempfile(pattern = paste0(".", basename(path), ".tmp-"), tmpdir = parent)
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  writeLines(lines, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) stop("Atomic rename failed for: ", path)
  invisible(path)
}

format_sci_math <- function(x, digits = 2L) {
  if (!is.finite(x)) return("--")
  if (x == 0) return("0")
  exponent <- floor(log10(abs(x)))
  mantissa <- x / 10^exponent
  sprintf(paste0("%.", digits, "f\\times 10^{%d}"), mantissa, exponent)
}

format_value_mcse <- function(value, se) {
  if (!is.finite(value)) return("--")
  if (!is.finite(se)) return(paste0("$", format_sci_math(value), "$"))
  paste0("$", format_sci_math(value), "\\;( ", format_sci_math(se), " )$")
}

format_pct_mcse <- function(value, se) {
  if (!is.finite(value)) return("--")
  if (!is.finite(se)) return(sprintf("%+.2f", value))
  sprintf("%+.2f (%.2f)", value, se)
}

order_to_plan <- function(x, source) {
  planned <- planned_grid()
  index <- match(config_key(planned), config_key(x))
  if (anyNA(index)) stop(source, " is missing a planned configuration")
  x[index, , drop = FALSE]
}

validate_table_inputs <- function(representative_abs, representative_pair, across, failures) {
  common <- c("kernel", "bw_multiplier")
  require_columns(
    representative_abs,
    c(common, "n", "N", "m", "reps_planned", "reps_ok",
      "mean_imse_sq", "mcse_imse_sq", "mean_imse", "mcse_imse"),
    "summary_by_cell data"
  )
  require_columns(
    representative_pair,
    c(common, "n", "N", "m", "pairs_planned", "pairs_ok",
      "paired_mean_pct_imse_sq", "paired_mcse_pct_imse_sq"),
    "paired_summary_by_cell data"
  )
  require_columns(across, c(common, "cells_planned", "cells_ok"),
                  "summary_across_cells data")
  require_columns(failures, c(common, "fits_planned", "fits_ok", "failures"),
                  "failure_summary_by_config data")

  planned_keys <- config_key(planned_grid())
  inputs <- list(
    summary_by_cell = representative_abs,
    paired_summary_by_cell = representative_pair,
    summary_across_cells = across,
    failure_summary_by_config = failures
  )
  for (name in names(inputs)) {
    x <- inputs[[name]]
    require_true(nrow(x) == 15L, paste0(name, " must have exactly 15 rows"))
    keys <- config_key(x)
    require_true(!anyDuplicated(keys), paste0(name, " contains duplicate configurations"))
    require_true(setequal(keys, planned_keys),
                 paste0(name, " does not contain the exact planned configuration grid"))
  }

  require_true(
    all(representative_abs$n == 100 & representative_abs$N == 5 &
          representative_abs$m == 200),
    "Absolute summaries must be restricted to n=100, N=5, m=200"
  )
  require_true(
    all(representative_pair$n == 100 & representative_pair$N == 5 &
          representative_pair$m == 200),
    "Paired summaries must be restricted to n=100, N=5, m=200"
  )
  require_true(
    all(representative_abs$reps_planned == 200 &
          representative_abs$reps_ok == 200),
    "The representative cell must have 200 successful fits per configuration"
  )
  require_true(
    all(representative_pair$pairs_planned == 200 &
          representative_pair$pairs_ok == 200),
    "The representative cell must have 200 successful configuration--default pairs"
  )
  require_true(
    all(across$cells_planned == 20 & integer_like(across$cells_ok) &
          across$cells_ok >= 0 & across$cells_ok <= 20),
    "Across-cell summaries must describe the 20-cell full design"
  )
  require_true(
    all(failures$fits_planned == 4000 & integer_like(failures$fits_ok) &
          integer_like(failures$failures) & failures$fits_ok >= 0 &
          failures$failures >= 0 &
          failures$fits_ok + failures$failures == failures$fits_planned),
    "Failure summaries must account for all 4,000 planned fits per configuration"
  )

  displayed <- c(
    representative_abs$mean_imse_sq, representative_abs$mcse_imse_sq,
    representative_abs$mean_imse, representative_abs$mcse_imse,
    representative_pair$paired_mean_pct_imse_sq,
    representative_pair$paired_mcse_pct_imse_sq
  )
  require_true(all(is.finite(displayed)), "A displayed estimate or MCSE is non-finite")
  invisible(TRUE)
}

write_sensitivity_table <- function(representative_abs, representative_pair,
                                    across, failures, path) {
  validate_table_inputs(representative_abs, representative_pair, across, failures)
  representative_abs <- order_to_plan(representative_abs, "summary_by_cell data")
  representative_pair <- order_to_plan(representative_pair, "paired_summary_by_cell data")
  across <- order_to_plan(across, "summary_across_cells data")
  failures <- order_to_plan(failures, "failure_summary_by_config data")
  planned <- planned_grid()

  labels <- c(gauss = "Gaussian", epan = "Epanechnikov", quar = "Quartic")
  body <- character(nrow(planned))
  for (i in seq_len(nrow(planned))) {
    is_default <- planned$kernel[i] == "gauss" &&
      abs(planned$bw_multiplier[i] - 1) < 1e-12
    fields <- c(
      unname(labels[planned$kernel[i]]),
      sprintf("%.2f", planned$bw_multiplier[i]),
      format_value_mcse(representative_abs$mean_imse_sq[i],
                        representative_abs$mcse_imse_sq[i]),
      format_pct_mcse(representative_pair$paired_mean_pct_imse_sq[i],
                      representative_pair$paired_mcse_pct_imse_sq[i]),
      format_value_mcse(representative_abs$mean_imse[i],
                        representative_abs$mcse_imse[i]),
      sprintf("%d/20", as.integer(across$cells_ok[i])),
      as.character(as.integer(failures$failures[i]))
    )
    if (is_default) fields <- paste0("\\textbf{", fields, "}")
    body[i] <- paste0(paste(fields, collapse = " & "), " \\\\")
  }

  lines <- c(
    "% Parentheses contain Monte Carlo standard errors.",
    "\\begin{table}[tbp]",
    "\\centering",
    "\\caption{Sensitivity of the proposed method to the temporal smoothing kernel and joint bandwidth multiplier $c_h$ under Setting~1.  Absolute results are shown for the prespecified representative cell $n=100$, $N=5$, and $m=200$.  The paired percentage change is the Monte Carlo mean of the replication-level changes relative to the default.  All configurations had 200 successful fits and 200 successful configuration--default pairs in this representative cell.  The valid-cells column gives the number of the 20 $(n,N)$ cells with at least one successful configuration--default pair.  Failures are counted over all 4,000 planned fits per configuration.  Parentheses contain Monte Carlo standard errors, and bold type identifies the untouched Gaussian default.}",
    "\\label{tab:sensitivity}",
    "\\small",
    "\\setlength{\\tabcolsep}{3.2pt}",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{llccccc}",
    "\\toprule",
    "Kernel & $c_h$ & IMSE (MCSE) & Paired change (\\%; MCSE) & $\\int d_W$ (MCSE) & Valid cells & Failures \\\\",
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\end{table}",
    ""
  )
  atomic_write_lines(lines, path)
}

read_summary <- function(path) {
  if (!file.exists(path)) stop("Required summary file is missing: ", path)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

validate_full_summaries <- function(summary_abs, summary_pair, across, failures) {
  require_true(nrow(summary_abs) == 300L,
               "summary_by_cell.csv must contain 300 rows")
  require_true(nrow(summary_pair) == 300L,
               "paired_summary_by_cell.csv must contain 300 rows")
  require_true(nrow(across) == 15L,
               "summary_across_cells.csv must contain 15 rows")
  require_true(nrow(failures) == 15L,
               "failure_summary_by_config.csv must contain 15 rows")

  for (x in list(summary_abs, summary_pair, across, failures)) {
    require_columns(x, c("run_id", "setting", "kernel", "bw_multiplier"),
                    "summary file")
    require_true(length(unique(x$run_id)) == 1L && nzchar(unique(x$run_id)),
                 "Each summary file must contain one nonempty run_id")
    require_true(all(x$setting == 1), "All summary rows must use Setting 1")
  }
  run_ids <- vapply(
    list(summary_abs, summary_pair, across, failures),
    function(x) unique(x$run_id), character(1L)
  )
  require_true(length(unique(run_ids)) == 1L,
               "Summary files do not have the same run_id")

  expected_cells <- expand.grid(
    n = c(20L, 50L, 100L, 200L, 400L),
    N = c(3L, 5L, 10L, 20L),
    m = 200L,
    KEEP.OUT.ATTRS = FALSE
  )
  cell_key <- function(x) paste(x$n, x$N, x$m, sep = "|")
  planned_keys <- config_key(planned_grid())
  for (item in list(summary_abs, summary_pair)) {
    require_columns(item, c("n", "N", "m"), "cell-level summary file")
    require_true(
      !anyDuplicated(paste(cell_key(item), config_key(item), sep = "|")),
      "A cell-level summary contains duplicate cell/configuration rows"
    )
    require_true(setequal(unique(cell_key(item)), cell_key(expected_cells)),
                 "A cell-level summary does not contain the exact 20-cell design")
    require_true(setequal(unique(config_key(item)), planned_keys),
                 "A cell-level summary does not contain the exact configuration grid")
  }
  invisible(TRUE)
}

generate_sensitivity_table <- function(summary_dir, table_out) {
  if (!dir.exists(summary_dir)) stop("Summary directory does not exist: ", summary_dir)
  summary_dir <- normalizePath(summary_dir, mustWork = TRUE)
  summary_abs <- read_summary(file.path(summary_dir, "summary_by_cell.csv"))
  summary_pair <- read_summary(file.path(summary_dir, "paired_summary_by_cell.csv"))
  across <- read_summary(file.path(summary_dir, "summary_across_cells.csv"))
  failures <- read_summary(file.path(summary_dir, "failure_summary_by_config.csv"))
  validate_full_summaries(summary_abs, summary_pair, across, failures)

  representative_abs <- summary_abs[
    summary_abs$n == 100 & summary_abs$N == 5 & summary_abs$m == 200,
    , drop = FALSE
  ]
  representative_pair <- summary_pair[
    summary_pair$n == 100 & summary_pair$N == 5 & summary_pair$m == 200,
    , drop = FALSE
  ]
  write_sensitivity_table(
    representative_abs, representative_pair, across, failures, table_out
  )
}

main <- function() {
  kv <- parse_kv_args(commandArgs(trailingOnly = TRUE))
  unknown <- setdiff(names(kv), c("summary_dir", "table_out"))
  if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
  if (is.null(kv$summary_dir) || !nzchar(kv$summary_dir)) {
    stop("Required argument: summary_dir=PATH")
  }
  if (is.null(kv$table_out) || !nzchar(kv$table_out)) {
    stop("Required argument: table_out=PATH")
  }
  generate_sensitivity_table(kv$summary_dir, kv$table_out)
  cat("LaTeX table:", normalizePath(kv$table_out, mustWork = TRUE), "\n")
}

if (sys.nframe() == 0L) main()
