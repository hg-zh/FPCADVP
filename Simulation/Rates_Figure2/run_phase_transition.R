#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressMessages(library(parallel))

# ---- RStudio run config ----
# Edit this block when running the script directly via Source / Run in RStudio.
# Command-line key=value arguments still override these values.
run_cfg <- list(
  project_dir = NULL,
  outdir = "outputs",
  label = "densegrid",
  n_vec = c(20L, 30L, 40L, 50L, 70L, 100L, 150L),
  N_vec = c(3L, 4L, 5L, 7L, 10L),
  m_vec = c(200L),
  reps = 200L,
  seed_base = 20260323L,
  cores = 5L,
  t_design = "random",
  latent = TRUE,
  methodXi = "CE",
  decayR = 1,
  eps_u = 0.05,
  u_cap = 0.95,
  Z_cap = 8,
  t_grid_len = 51L,
  slope_n_min = 70L,
  save_each_rep = FALSE,
  mod_script = NULL,
  raw_file = "raw_results_densegrid.rds",
  summary_file = "summary_densegrid.csv",
  diagnostics_file = "diagnostics_densegrid.csv",
  session_file = "sessionInfo.txt",
  notes_file = "notes.md",
  config_file = "run_config_densegrid.rds"
)

parse_kv_args <- function(args) {
  out <- list()
  if (length(args) == 0) return(out)
  for (a in args) {
    if (!grepl("=", a, fixed = TRUE)) next
    key <- sub("=.*$", "", a)
    val <- sub("^[^=]*=", "", a)
    out[[key]] <- val
  }
  out
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0) return(default)
  if (is.logical(x)) return(isTRUE(x))
  s <- tolower(trimws(as.character(x)[1]))
  if (s %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (s %in% c("false", "f", "0", "no", "n")) return(FALSE)
  default
}

as_int_or_default <- function(x, default) {
  if (is.null(x)) return(as.integer(default))
  out <- suppressWarnings(as.integer(x))
  if (!is.finite(out)) return(as.integer(default))
  out
}

as_num_or_default <- function(x, default) {
  if (is.null(x)) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(x))
  if (!is.finite(out)) return(as.numeric(default))
  out
}

parse_int_vec <- function(x) {
  x <- paste(as.character(x), collapse = ",")
  x <- gsub("\\s+", "", x)
  if (nchar(x) == 0) return(integer(0))
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

as_method_or_null <- function(x) {
  if (is.null(x)) return(NULL)
  s <- trimws(as.character(x)[1])
  if (tolower(s) %in% c("", "auto", "null", "none")) return(NULL)
  s
}

fmt_num <- function(x, digits = 4) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
}

format_seconds <- function(x) {
  x <- as.numeric(x)
  if (!is.finite(x)) return("NA")
  sprintf("%.1f seconds (%.2f minutes)", x, x / 60)
}

call_with_supported_args <- function(fn, args) {
  keep <- names(args) %in% names(formals(fn))
  do.call(fn, args[keep])
}

resolve_outfile <- function(outdir, file_arg, default_basename) {
  path <- file_arg %||% default_basename
  if (grepl("^/", path)) return(path)
  file.path(outdir, path)
}

detect_mod_script <- function(user_path = NULL) {
  if (!is.null(user_path)) return(user_path)
  candidates <- c(
    file.path( "transport_fpca_modular.R"),
    file.path( "transport_fpca_modular.R")
  )
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0) {
    stop("Could not find transport_fpca_modular.R under proposed_method/ or codes/.")
  }
  hits[1]
}

flatten_numeric <- function(x) {
  if (is.null(x)) return(numeric(0))
  if (is.list(x)) {
    return(unlist(lapply(x, function(y) as.numeric(as.matrix(y))), use.names = FALSE))
  }
  as.numeric(as.matrix(x))
}

transport_inner_boundary_frac <- function(T_obj, tol = 1e-6) {
  mats <- if (is.list(T_obj)) lapply(T_obj, as.matrix) else list(as.matrix(T_obj))
  vals <- unlist(lapply(mats, function(M) {
    if (ncol(M) <= 2) return(numeric(0))
    as.numeric(M[, 2:(ncol(M) - 1), drop = FALSE])
  }), use.names = FALSE)
  if (length(vals) == 0) return(NA_real_)
  mean(vals <= tol | vals >= 1 - tol)
}

frac_at_cap <- function(x, cap, tol = 1e-8) {
  if (!is.finite(cap)) return(NA_real_)
  vals <- flatten_numeric(x)
  if (length(vals) == 0) return(NA_real_)
  mean(abs(vals) >= cap - tol)
}

range_or_na <- function(x) {
  vals <- flatten_numeric(x)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(c(NA_real_, NA_real_))
  range(vals)
}

safe_var <- function(x) {
  x <- flatten_numeric(x)
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::var(x)
}

interp_Tmat_time <- function(Th, t_from, t_to) {
  Th <- as.matrix(Th)
  out <- matrix(NA_real_, nrow = length(t_to), ncol = ncol(Th))
  for (j in seq_len(ncol(Th))) {
    out[, j] <- approx(x = t_from, y = Th[, j], xout = t_to, rule = 2)$y
  }
  out
}

summarize_cells <- function(ok_df) {
  split_key <- interaction(ok_df$n, ok_df$N, ok_df$m, drop = TRUE)
  groups <- split(ok_df, split_key)
  diag_cols <- setdiff(names(ok_df), c("rep", "n", "N", "m", "iw1", "imse", "ok", "err"))
  out <- lapply(groups, function(d) {
    diag_means <- if (length(diag_cols) > 0) {
      stats::setNames(
        as.list(vapply(diag_cols, function(col) mean(d[[col]], na.rm = TRUE), numeric(1))),
        paste0("mean_", diag_cols)
      )
    } else {
      list()
    }
    data.frame(
      n = d$n[1],
      N = d$N[1],
      m = d$m[1],
      reps_ok = nrow(d),
      mean_imse = mean(d$imse),
      sd_imse = if (nrow(d) > 1) stats::sd(d$imse) else 0,
      se_imse = if (nrow(d) > 1) stats::sd(d$imse) / sqrt(nrow(d)) else 0,
      mean_iw1 = mean(d$iw1),
      sd_iw1 = if (nrow(d) > 1) stats::sd(d$iw1) else 0,
      se_iw1 = if (nrow(d) > 1) stats::sd(d$iw1) / sqrt(nrow(d)) else 0,
      stringsAsFactors = FALSE
    ) |>
      cbind(as.data.frame(diag_means, stringsAsFactors = FALSE))
  })
  out <- do.call(rbind, out)
  out[order(out$N, out$n, out$m), , drop = FALSE]
}

estimate_slopes_raw <- function(ok_df, slope_n_min, theoretical_slope, a_poly) {
  slope_df <- ok_df[ok_df$n >= slope_n_min & ok_df$imse > 0 & is.finite(ok_df$imse), , drop = FALSE]
  if (nrow(slope_df) == 0) return(data.frame())
  split_key <- interaction(slope_df$N, slope_df$m, drop = TRUE)
  groups <- split(slope_df, split_key)
  out <- lapply(groups, function(d) {
    fit <- stats::lm(log(imse) ~ log(n), data = d)
    coef_tab <- summary(fit)$coefficients
    data.frame(
      N = d$N[1],
      m = d$m[1],
      slope_tail_n_values = paste(sort(unique(d$n)), collapse = ","),
      slope_fit_n_min = slope_n_min,
      slope_est = coef_tab["log(n)", "Estimate"],
      slope_se = coef_tab["log(n)", "Std. Error"],
      theoretical_slope = theoretical_slope,
      a_poly = a_poly,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

estimate_curve_slopes <- function(cell_df, slope_n_min) {
  if (nrow(cell_df) == 0) return(data.frame())
  split_key <- interaction(cell_df$N, cell_df$m, drop = TRUE)
  groups <- split(cell_df, split_key)
  out <- lapply(groups, function(d) {
    fit_all <- stats::lm(log(mean_imse) ~ log(n), data = d)
    d_tail <- d[d$n >= slope_n_min, , drop = FALSE]
    fit_tail <- if (nrow(d_tail) >= 2) stats::lm(log(mean_imse) ~ log(n), data = d_tail) else NULL
    data.frame(
      N = d$N[1],
      m = d$m[1],
      curve_slope_all = unname(stats::coef(fit_all)[["log(n)"]]),
      curve_slope_tail = if (!is.null(fit_tail)) unname(stats::coef(fit_tail)[["log(n)"]]) else NA_real_,
      curve_tail_n_values = if (nrow(d_tail) > 0) paste(sort(unique(d_tail$n)), collapse = ",") else "",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

diagnostics_subset <- function(cell_df) {
  keep <- c(
    "n", "N", "m", "reps_ok",
    "mean_sample_var", "mean_u_obs_var",
    "mean_true_min", "mean_true_max",
    "mean_obs_min", "mean_obs_max",
    "mean_pred_min", "mean_pred_max",
    "mean_frac_true_inner_boundary",
    "mean_frac_obs_inner_boundary",
    "mean_frac_pred_inner_boundary",
    "mean_max_abs_u_raw",
    "mean_frac_abs_u_raw_gt1",
    "mean_max_abs_u_used",
    "mean_frac_u_obs_at_cap",
    "mean_frac_u_grid_at_cap",
    "mean_frac_z_grid_at_cap",
    "mean_xgrid_len",
    "mean_eval_grid_len",
    "mean_workgrid_len",
    "mean_kappa_true"
  )
  keep <- keep[keep %in% names(cell_df)]
  cell_df[, keep, drop = FALSE]
}

write_notes <- function(notes_path, config, cell_df, raw_slope_df, curve_slope_df, runtime_sec, extra_lines = character(0)) {
  err_vec <- config$error_messages
  err_lines <- if (length(err_vec) > 0) {
    paste(sprintf("- `%s`: %d", names(err_vec), as.integer(err_vec)), collapse = "\n")
  } else {
    "- No failed fits."
  }

  raw_slope_lines <- if (nrow(raw_slope_df) > 0) {
    paste(apply(raw_slope_df[order(raw_slope_df$N), , drop = FALSE], 1, function(row) {
      sprintf("| %s | %s | %s | %s |",
              row[["N"]],
              row[["slope_tail_n_values"]],
              fmt_num(as.numeric(row[["slope_est"]]), digits = 4),
              fmt_num(as.numeric(row[["slope_se"]]), digits = 4))
    }), collapse = "\n")
  } else {
    "| NA | NA | NA | NA |"
  }

  curve_slope_lines <- if (nrow(curve_slope_df) > 0) {
    paste(apply(curve_slope_df[order(curve_slope_df$N), , drop = FALSE], 1, function(row) {
      sprintf("| %s | %s | %s |",
              row[["N"]],
              fmt_num(as.numeric(row[["curve_slope_all"]]), digits = 4),
              fmt_num(as.numeric(row[["curve_slope_tail"]]), digits = 4))
    }), collapse = "\n")
  } else {
    "| NA | NA | NA |"
  }

  diag_lines <- character(0)
  if (nrow(cell_df) > 0) {
    diag_lines <- c(
      sprintf("- IMSE was evaluated on the fixed grid `t_grid_out = seq(0, 1, length.out = %d)` for every setting.", config$t_grid_len),
      sprintf("- Mean generated sample variance across cells ranged from `%s` to `%s`.", fmt_num(min(cell_df$mean_sample_var, na.rm = TRUE), 4), fmt_num(max(cell_df$mean_sample_var, na.rm = TRUE), 4)),
      sprintf("- True transport values stayed in [`%s`, `%s`]; predicted transport values stayed in [`%s`, `%s`].",
              fmt_num(min(cell_df$mean_true_min, na.rm = TRUE), 4),
              fmt_num(max(cell_df$mean_true_max, na.rm = TRUE), 4),
              fmt_num(min(cell_df$mean_pred_min, na.rm = TRUE), 4),
              fmt_num(max(cell_df$mean_pred_max, na.rm = TRUE), 4)),
      sprintf("- Mean clamp / saturation diagnostics: `frac_abs_u_raw_gt1` in [`%s`, `%s`], `frac_u_grid_at_cap` in [`%s`, `%s`], `frac_pred_inner_boundary` in [`%s`, `%s`].",
              fmt_num(min(cell_df$mean_frac_abs_u_raw_gt1, na.rm = TRUE), 4),
              fmt_num(max(cell_df$mean_frac_abs_u_raw_gt1, na.rm = TRUE), 4),
              fmt_num(min(cell_df$mean_frac_u_grid_at_cap, na.rm = TRUE), 4),
              fmt_num(max(cell_df$mean_frac_u_grid_at_cap, na.rm = TRUE), 4),
              fmt_num(min(cell_df$mean_frac_pred_inner_boundary, na.rm = TRUE), 4),
              fmt_num(max(cell_df$mean_frac_pred_inner_boundary, na.rm = TRUE), 4))
    )
  }

  lines <- c(
    "# Phase Transition Run Notes",
    "",
    "## Theory Mapping",
    sprintf("- Corollary 1 uses slope `-(a-1)/(3a+1)` with polynomial eigenvalue decay `lambda_k ~ k^{-a}`."),
    sprintf("- Here `decayR = %.2f` means score SDs decay like `k^{-decayR}`, so variances decay like `k^{-2*decayR}` and `a = %.1f`.", config$decayR, config$a_poly),
    sprintf("- The theoretical log-log slope is `%s`.", fmt_num(config$theoretical_slope, digits = 4)),
    "",
    "## Final Settings",
    sprintf("- Grid: `n = c(%s)`, `N = c(%s)`, `m = c(%s)`.", paste(config$n_vec, collapse = ", "), paste(config$N_vec, collapse = ", "), paste(config$m_vec, collapse = ", ")),
    sprintf("- Replications: `%d`; seed base: `%d`; cores: `%d`.", config$reps, config$seed_base, config$cores),
    sprintf("- Design: `%s`; latent = `%s`; methodXi = `%s`.", config$t_design, config$latent, config$methodXi),
    sprintf("- Stability knobs: `eps_u = %.2f`, `u_cap = %.2f`, `Z_cap = %.1f`, `t_grid_len = %d`.", config$eps_u, config$u_cap, config$Z_cap, config$t_grid_len),
    sprintf("- Runtime: %s.", format_seconds(runtime_sec)),
    sprintf("- Modular method file sourced from `%s`.", config$mod_script),
    "",
    "## Tail Slopes From Replication-Level IMSE",
    sprintf("- Slopes use `n >= %d`.", config$slope_n_min),
    "",
    "| N | tail n values | slope | SE |",
    "| --- | --- | --- | --- |",
    raw_slope_lines,
    "",
    "## Mean-Curve Slopes For Plotting",
    "| N | all-n curve slope | tail-only curve slope |",
    "| --- | --- | --- |",
    curve_slope_lines,
    "",
    "## Diagnostics",
    diag_lines,
    "",
    "## Flattening Diagnosis",
    extra_lines,
    "",
    "## Failure Check",
    err_lines
  )

  writeLines(lines, con = notes_path)
}

kv_cli <- parse_kv_args(commandArgs(trailingOnly = TRUE))
kv <- if (length(kv_cli) == 0) run_cfg else utils::modifyList(run_cfg, kv_cli)

if (!is.null(kv$project_dir) && nzchar(as.character(kv$project_dir))) {
  setwd(as.character(kv$project_dir))
}

default_n_vec <- c(20L, 30L, 40L, 50L, 70L, 100L, 150L, 200L)
default_N_vec <- c(2L, 3L, 4L, 5L, 7L, 10L)
default_m_vec <- c(200L)

n_vec <- if (!is.null(kv$n_vec)) parse_int_vec(kv$n_vec) else default_n_vec
N_vec <- if (!is.null(kv$N_vec)) parse_int_vec(kv$N_vec) else default_N_vec
m_vec <- if (!is.null(kv$m_vec)) parse_int_vec(kv$m_vec) else default_m_vec

if (length(n_vec) == 0 || any(!is.finite(n_vec))) stop("n_vec must be a non-empty integer vector.")
if (length(N_vec) == 0 || any(!is.finite(N_vec))) stop("N_vec must be a non-empty integer vector.")
if (length(m_vec) == 0 || any(!is.finite(m_vec))) stop("m_vec must be a non-empty integer vector.")

reps <- as_int_or_default(kv$reps, 120L)
seed_base <- as_int_or_default(kv$seed_base, 20260320L)
cores_user <- as_int_or_default(kv$cores, min(5L, parallel::detectCores(logical = TRUE)))
cores <- max(1L, min(5L, cores_user))
outdir <- kv$outdir %||% "outputs"
label <- kv$label %||% "densegrid"
mod_script <- detect_mod_script(kv$mod_script)
t_design <- kv$t_design %||% "random"
latent <- as_bool(kv$latent, TRUE)
methodXi <- as_method_or_null(kv$methodXi %||% "CE")
decayR <- as_num_or_default(kv$decayR, 1)
eps_u <- as_num_or_default(kv$eps_u, 0.05)
u_cap <- as_num_or_default(kv$u_cap, 0.95)
Z_cap <- as_num_or_default(kv$Z_cap, 8)
t_grid_len <- as_int_or_default(kv$t_grid_len, 51L)
slope_n_min <- as_int_or_default(kv$slope_n_min, 70L)
save_each_rep <- as_bool(kv$save_each_rep, FALSE)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

raw_path <- resolve_outfile(outdir, kv$raw_file, "raw_results_densegrid.rds")
summary_path <- resolve_outfile(outdir, kv$summary_file, "summary_densegrid.csv")
diag_path <- resolve_outfile(outdir, kv$diagnostics_file, "diagnostics_densegrid.csv")
session_path <- resolve_outfile(outdir, kv$session_file, "sessionInfo.txt")
notes_path <- resolve_outfile(outdir, kv$notes_file, "notes.md")

if (!t_design %in% c("random", "fixed")) stop("t_design must be 'random' or 'fixed'.")
if (!is.finite(u_cap) || u_cap <= 0 || u_cap >= 1) stop("u_cap must be in (0, 1).")
if (!is.finite(eps_u) || eps_u <= 0 || eps_u >= 1) stop("eps_u must be in (0, 1).")
if (!is.finite(Z_cap) || Z_cap <= 0) stop("Z_cap must be positive.")
if (!is.finite(decayR) || decayR <= 0) stop("decayR must be positive.")
if (reps < 1) stop("reps must be >= 1.")
if (t_grid_len < 2) stop("t_grid_len must be >= 2.")

t_grid_out <- seq(0, 1, length.out = t_grid_len)
a_poly <- 2 * decayR
theoretical_slope <- -(a_poly - 1) / (3 * a_poly + 1)

check_env <- new.env(parent = globalenv())
sys.source(mod_script, envir = check_env)
required_funs <- c("simulate_transport_process", "fit_transport_fpca", "true_transports_on_grid", "trapz_rows", "trapz_vec")
missing_funs <- required_funs[!vapply(required_funs, exists, logical(1), envir = check_env, inherits = FALSE)]
if (length(missing_funs) > 0) {
  stop("mod_script is missing required functions: ", paste(missing_funs, collapse = ", "))
}
rm(check_env)

run_one_rep <- function(rep_id,
                        n_vec, N_vec, m_vec,
                        seed_base,
                        t_design,
                        latent, methodXi,
                        eps_u, u_cap, Z_cap,
                        t_grid_out,
                        decayR,
                        outdir,
                        save_each_rep) {
  rows <- vector("list", length(n_vec) * length(N_vec) * length(m_vec))
  idx <- 0L

  for (n in n_vec) {
    for (N in N_vec) {
      for (m in m_vec) {
        seed <- seed_base + rep_id * 1000000L + n * 1000L + N * 10L + m

        out <- tryCatch({
          sim <- simulate_transport_process(
            n = n,
            N = N,
            m = m,
            t_design = t_design,
            seed = seed,
            decayR = decayR
          )

          fit_args <- list(
            T_list = sim$T_hat,
            t_list = sim$t_list,
            xgrid = sim$xgrid,
            kappa = sim$kappa_true,
            latent = latent,
            methodXi = methodXi,
            fpca_optns = list(methodMuCovEst = "smooth"),
            eps_u = eps_u,
            u_cap = u_cap,
            Z_cap = Z_cap,
            project = TRUE,
            verbose = FALSE
          )
          fit <- call_with_supported_args(fit_transport_fpca, fit_args)

          T_true <- lapply(seq_len(n), function(i) {
            true_transports_on_grid(
              Xi_row = sim$Xi[i, ],
              T0_row = sim$T0_true[i, ],
              t_grid = t_grid_out,
              xgrid = sim$xgrid,
              K = sim$K,
              g = sim$g,
              project = TRUE
            )
          })

          T_pred <- fit$T_hat
          if (length(fit$workGrid) != length(t_grid_out) ||
              max(abs(fit$workGrid - t_grid_out)) > 1e-12) {
            T_pred <- lapply(T_pred, function(Th) interp_Tmat_time(Th, fit$workGrid, t_grid_out))
          }

          iw1_i <- numeric(n)
          imse_i <- numeric(n)
          for (i in seq_len(n)) {
            Th <- as.matrix(T_pred[[i]])
            Tt <- as.matrix(T_true[[i]])
            if (!all(dim(Th) == dim(Tt))) stop("dim mismatch after interpolation")
            w1_t <- trapz_rows(sim$xgrid, abs(Th - Tt))
            iw1_i[i] <- trapz_vec(t_grid_out, w1_t)
            imse_i[i] <- trapz_vec(t_grid_out, w1_t^2)
          }

          sample_vals <- flatten_numeric(sim$X_samp)
          true_rng <- range_or_na(sim$T_true)
          obs_rng <- range_or_na(sim$T_hat)
          pred_rng <- range_or_na(T_pred)

          list(
            ok = TRUE,
            iw1 = mean(iw1_i),
            imse = mean(imse_i),
            sample_var = safe_var(sample_vals),
            u_obs_var = safe_var(sim$U_obs),
            true_min = true_rng[1],
            true_max = true_rng[2],
            obs_min = obs_rng[1],
            obs_max = obs_rng[2],
            pred_min = pred_rng[1],
            pred_max = pred_rng[2],
            frac_true_inner_boundary = transport_inner_boundary_frac(sim$T_true),
            frac_obs_inner_boundary = transport_inner_boundary_frac(sim$T_hat),
            frac_pred_inner_boundary = transport_inner_boundary_frac(T_pred),
            max_abs_u_raw = fit$diagnostics$max_abs_U_raw,
            frac_abs_u_raw_gt1 = fit$diagnostics$frac_abs_U_raw_gt1,
            max_abs_u_used = fit$diagnostics$max_abs_U_used,
            frac_u_obs_at_cap = frac_at_cap(fit$U_hat_obs, fit$u_cap),
            frac_u_grid_at_cap = frac_at_cap(fit$U_hat, fit$u_cap),
            frac_z_grid_at_cap = if (is.matrix(fit$Z_hat) && is.finite(Z_cap)) mean(abs(fit$Z_hat) >= Z_cap - 1e-8) else NA_real_,
            xgrid_len = length(sim$xgrid),
            eval_grid_len = length(t_grid_out),
            workgrid_len = length(fit$workGrid),
            kappa_true = sim$kappa_true
          )
        }, error = function(e) {
          list(
            ok = FALSE,
            iw1 = NA_real_,
            imse = NA_real_,
            sample_var = NA_real_,
            u_obs_var = NA_real_,
            true_min = NA_real_,
            true_max = NA_real_,
            obs_min = NA_real_,
            obs_max = NA_real_,
            pred_min = NA_real_,
            pred_max = NA_real_,
            frac_true_inner_boundary = NA_real_,
            frac_obs_inner_boundary = NA_real_,
            frac_pred_inner_boundary = NA_real_,
            max_abs_u_raw = NA_real_,
            frac_abs_u_raw_gt1 = NA_real_,
            max_abs_u_used = NA_real_,
            frac_u_obs_at_cap = NA_real_,
            frac_u_grid_at_cap = NA_real_,
            frac_z_grid_at_cap = NA_real_,
            xgrid_len = NA_real_,
            eval_grid_len = NA_real_,
            workgrid_len = NA_real_,
            kappa_true = NA_real_,
            err = conditionMessage(e)
          )
        })

        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          rep = rep_id,
          n = n,
          N = N,
          m = m,
          iw1 = out$iw1,
          imse = out$imse,
          sample_var = out$sample_var,
          u_obs_var = out$u_obs_var,
          true_min = out$true_min,
          true_max = out$true_max,
          obs_min = out$obs_min,
          obs_max = out$obs_max,
          pred_min = out$pred_min,
          pred_max = out$pred_max,
          frac_true_inner_boundary = out$frac_true_inner_boundary,
          frac_obs_inner_boundary = out$frac_obs_inner_boundary,
          frac_pred_inner_boundary = out$frac_pred_inner_boundary,
          max_abs_u_raw = out$max_abs_u_raw,
          frac_abs_u_raw_gt1 = out$frac_abs_u_raw_gt1,
          max_abs_u_used = out$max_abs_u_used,
          frac_u_obs_at_cap = out$frac_u_obs_at_cap,
          frac_u_grid_at_cap = out$frac_u_grid_at_cap,
          frac_z_grid_at_cap = out$frac_z_grid_at_cap,
          xgrid_len = out$xgrid_len,
          eval_grid_len = out$eval_grid_len,
          workgrid_len = out$workgrid_len,
          kappa_true = out$kappa_true,
          ok = out$ok,
          err = if (!isTRUE(out$ok)) out$err else "",
          stringsAsFactors = FALSE
        )

        gc(verbose = FALSE)
      }
    }
  }

  df <- do.call(rbind, rows)
  if (isTRUE(save_each_rep)) {
    saveRDS(df, file.path(outdir, sprintf("rep_%03d.rds", rep_id)))
  }
  df
}

start_time <- Sys.time()
cat(sprintf("Starting phase-transition run at %s\n", format(start_time, tz = Sys.timezone(), usetz = TRUE)))
cat(sprintf("label=%s | reps=%d | cores=%d | outdir=%s\n", label, reps, cores, outdir))
cat(sprintf("grid: n=%s | N=%s | m=%s\n", paste(n_vec, collapse = ","), paste(N_vec, collapse = ","), paste(m_vec, collapse = ",")))

cl <- makeCluster(cores)
on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)

clusterExport(cl, varlist = c("mod_script"), envir = environment())
invisible(clusterEvalQ(cl, {
  suppressMessages(library(fdapace))
  source(mod_script)
}))

clusterExport(
  cl,
  varlist = c("call_with_supported_args", "interp_Tmat_time", "run_one_rep",
              "flatten_numeric", "transport_inner_boundary_frac", "frac_at_cap", "range_or_na", "safe_var",
              "n_vec", "N_vec", "m_vec", "seed_base", "t_design", "latent", "methodXi",
              "eps_u", "u_cap", "Z_cap", "t_grid_out", "decayR", "outdir", "save_each_rep"),
  envir = environment()
)

ptm <- proc.time()
res_list <- parLapply(cl, seq_len(reps), function(rep_id) {
  run_one_rep(
    rep_id = rep_id,
    n_vec = n_vec,
    N_vec = N_vec,
    m_vec = m_vec,
    seed_base = seed_base,
    t_design = t_design,
    latent = latent,
    methodXi = methodXi,
    eps_u = eps_u,
    u_cap = u_cap,
    Z_cap = Z_cap,
    t_grid_out = t_grid_out,
    decayR = decayR,
    outdir = outdir,
    save_each_rep = save_each_rep
  )
})
run_time <- proc.time() - ptm
stopCluster(cl)

raw_df <- do.call(rbind, res_list)
saveRDS(raw_df, raw_path)

ok_df <- raw_df[raw_df$ok %in% c(TRUE, 1), , drop = FALSE]
cell_df <- data.frame()
raw_slope_df <- data.frame()
curve_slope_df <- data.frame()

if (nrow(ok_df) == 0) {
  write.csv(data.frame(), summary_path, row.names = FALSE)
  write.csv(data.frame(), diag_path, row.names = FALSE)
} else {
  cell_df <- summarize_cells(ok_df)
  raw_slope_df <- estimate_slopes_raw(ok_df, slope_n_min = slope_n_min, theoretical_slope = theoretical_slope, a_poly = a_poly)
  curve_slope_df <- estimate_curve_slopes(cell_df, slope_n_min = slope_n_min)
  if (nrow(raw_slope_df) > 0) {
    cell_df <- merge(cell_df, raw_slope_df, by = c("N", "m"), all.x = TRUE, sort = FALSE)
  }
  if (nrow(curve_slope_df) > 0) {
    cell_df <- merge(cell_df, curve_slope_df, by = c("N", "m"), all.x = TRUE, sort = FALSE)
  }
  cell_df <- cell_df[order(cell_df$N, cell_df$n, cell_df$m), , drop = FALSE]
  write.csv(cell_df, summary_path, row.names = FALSE)
  write.csv(diagnostics_subset(cell_df), diag_path, row.names = FALSE)
}

err_vec <- raw_df$err[nzchar(raw_df$err)]
err_tab <- if (length(err_vec) > 0) sort(table(err_vec), decreasing = TRUE) else integer(0)

config <- list(
  label = label,
  reps = reps,
  seed_base = seed_base,
  cores = cores,
  t_design = t_design,
  latent = latent,
  methodXi = if (is.null(methodXi)) "AUTO" else methodXi,
  decayR = decayR,
  a_poly = a_poly,
  theoretical_slope = theoretical_slope,
  eps_u = eps_u,
  u_cap = u_cap,
  Z_cap = Z_cap,
  t_grid_len = t_grid_len,
  slope_n_min = slope_n_min,
  n_vec = n_vec,
  N_vec = N_vec,
  m_vec = m_vec,
  mod_script = normalizePath(mod_script),
  error_messages = err_tab
)
saveRDS(config, resolve_outfile(outdir, kv$config_file, paste0("run_config_", label, ".rds")))

note_lines <- c(
  "- This run writes versioned outputs so the earlier sparse-grid files remain available for comparison.",
  "- The main dense-grid plot focuses on `N <= 10` so that small-N curvature is visible and the earlier `N = 20` flattening does not dominate the figure."
)

write_notes(
  notes_path = notes_path,
  config = config,
  cell_df = cell_df,
  raw_slope_df = raw_slope_df,
  curve_slope_df = curve_slope_df,
  runtime_sec = as.numeric(run_time["elapsed"]),
  extra_lines = note_lines
)

session_lines <- capture.output(sessionInfo())
writeLines(session_lines, con = session_path)

end_time <- Sys.time()
cat(sprintf("Finished at %s\n", format(end_time, tz = Sys.timezone(), usetz = TRUE)))
cat(sprintf("Elapsed: %s\n", format_seconds(run_time["elapsed"])))
cat(sprintf("Saved raw results to %s\n", raw_path))
cat(sprintf("Saved summary to %s\n", summary_path))
cat(sprintf("Saved diagnostics to %s\n", diag_path))
