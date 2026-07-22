#!/usr/bin/env Rscript

## para_R2_baselines_mse_t.R
## ------------------------------------------------------------
## Parallel Monte Carlo simulation driver for two baseline methods:
##   (A) LQD + time-FPCA
##   (B) Tangent/log-map Wasserstein + time-FPCA
##
## EXTENSION vs para_R2_baselines.R:
##   - optionally compute and SAVE pointwise MSE(t) curves for each method
##     where MSE(t) = mean_i d_{W,1}^2(That_i(t), Ttrue_i(t)).
##   - robust boolean parsing (supports T/F)
##   - prefers SLURM_CPUS_PER_TASK (shared-memory) for cores.
##
## Outputs (in outdir):
##   - raw_results.rds
##   - summary.csv
##   - sessionInfo.txt
##   - (optional) MSE curves under outdir/mse_curves/{LQD,Tangent}/
##       MSE_curve_rep###_n##_N##_m##.csv
##       (optional) ... .pdf if plot_mse_t=TRUE
##
## Usage example:
##   module load r/4.5.0
##   Rscript para_R2_baselines_mse_t.R outdir=Res_R2_baselines_mse_t \
##     save_mse_t=TRUE plot_mse_t=FALSE t_grid_len=51 methodXi=CE
## ------------------------------------------------------------

options(stringsAsFactors = FALSE)

suppressMessages(library(parallel))

# ---- simple key=value parser ----
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
  if (s %in% c("true","t","1","yes","y")) return(TRUE)
  if (s %in% c("false","f","0","no","n")) return(FALSE)
  default
}

as_int1 <- function(x, default = NA_integer_) {
  if (is.null(x) || length(x) == 0) return(default)
  s <- trimws(as.character(x)[1])
  if (tolower(s) %in% c("", "na", "nan", "null", "none")) return(NA_integer_)
  suppressWarnings(as.integer(s))
}

parse_int_vec <- function(x) {
  x <- gsub("\\s+", "", as.character(x))
  if (nchar(x) == 0) return(integer(0))
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

kv <- parse_kv_args(commandArgs(trailingOnly = TRUE))

# ---- config (defaults) ----
reps       <- if (!is.null(kv$reps)) as.integer(kv$reps) else 200
seed_base  <- if (!is.null(kv$seed_base)) as.integer(kv$seed_base) else 2000
outdir     <- kv$outdir %||% "Res_R2_baselines"
t_design   <- kv$t_design %||% "random"   # "random" or "fixed"

# Baseline knobs
Ku         <- if (!is.null(kv$Ku)) as.integer(kv$Ku) else 10
basis_type <- kv$basis_type %||% "bspline" # "bspline" or "cosine"
bs_degree  <- if (!is.null(kv$bs_degree)) as.integer(kv$bs_degree) else 3

smooth_deriv <- if (!is.null(kv$smooth_deriv)) as_bool(kv$smooth_deriv, FALSE) else FALSE
spar         <- if (!is.null(kv$spar)) as.numeric(kv$spar) else NULL

# Output/eval grid
t_grid_len <- if (!is.null(kv$t_grid_len)) as.integer(kv$t_grid_len) else 51
methodXi   <- kv$methodXi %||% "CE"   # PACE/CE recommended
decayR     <- if (!is.null(kv$decayR)) as.numeric(kv$decayR) else 1
if (!is.finite(decayR) || decayR <= 0) decayR <- 1

save_each_rep <- if (!is.null(kv$save_each_rep)) as_bool(kv$save_each_rep, TRUE) else TRUE

# ---- NEW: MSE(t) curves ----
save_mse_t   <- if (!is.null(kv$save_mse_t)) as_bool(kv$save_mse_t, FALSE) else FALSE
plot_mse_t   <- if (!is.null(kv$plot_mse_t)) as_bool(kv$plot_mse_t, FALSE) else FALSE
mse_rep_only <- as_int1(kv$mse_rep_only, default = NA_integer_)

# Parameter grid (defaults match your current R2 grid)
# Override via n_vec=20,50,... etc if desired.
n_vec <- if (!is.null(kv$n_vec)) parse_int_vec(kv$n_vec) else c(20, 50, 100, 200, 400)
N_vec <- if (!is.null(kv$N_vec)) parse_int_vec(kv$N_vec) else c(3, 5, 10, 20)
m_vec <- if (!is.null(kv$m_vec)) parse_int_vec(kv$m_vec) else c(200)
beta_a <- if (!is.null(kv$beta_a)) parse_int_vec(kv$beta_a) else c(3,4)
beta_b <- if (!is.null(kv$beta_b)) parse_int_vec(kv$beta_b) else c(1,2)

if (length(n_vec) < 1 || any(!is.finite(n_vec))) stop("n_vec is empty or invalid")
if (length(N_vec) < 1 || any(!is.finite(N_vec))) stop("N_vec is empty or invalid")
if (length(m_vec) < 1 || any(!is.finite(m_vec))) stop("m_vec is empty or invalid")

# Output time grid for all methods
if (!is.finite(t_grid_len) || t_grid_len < 2) stop("t_grid_len must be >= 2")
t_grid_out <- seq(0, 1, length.out = t_grid_len)

# ---- cores: prefer CPUS_PER_TASK for shared-memory jobs ----
cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
if (!is.finite(cores) || cores < 1) cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_NTASKS")))
if ((!is.finite(cores) || cores < 1) && !is.null(kv$cores)) cores <- as.integer(kv$cores)
if (!is.finite(cores) || cores < 1) cores <- 20

cores_env <- Sys.getenv("SLURM_NTASKS")
cpus_env  <- Sys.getenv("SLURM_CPUS_PER_TASK")

# ---- locate scripts ----
mod_script  <- kv$mod_script  %||% "transport_fpca_modular.R"
base_script <- kv$base_script %||% "baseline_fpca_methods.R"

if (!file.exists(mod_script)) {
  stop(sprintf("Cannot find '%s'. Put transport_fpca_modular.R in the same directory or pass mod_script=PATH.", mod_script))
}
if (!file.exists(base_script)) {
  stop(sprintf("Cannot find '%s'. Put baseline_fpca_methods.R in the same directory or pass base_script=PATH.", base_script))
}

mod_script  <- normalizePath(mod_script)
base_script <- normalizePath(base_script)

# ---- preflight: detect whether simulate_transport_process supports decayR ----
sim_supports_decayR <- FALSE
{
  tmp_env <- new.env(parent = globalenv())
  try({
    sys.source(mod_script, envir = tmp_env)
    if (exists("simulate_transport_process", envir = tmp_env, inherits = FALSE)) {
      ff <- get("simulate_transport_process", envir = tmp_env)
      sim_supports_decayR <- ("decayR" %in% names(formals(ff)))
    }
  }, silent = TRUE)
  rm(tmp_env)
}

# ---- output dir ----
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Curve output dirs
curve_root <- kv$curve_root %||% file.path(outdir, "mse_curves")
curve_dir_lqd <- file.path(curve_root, "LQD")
curve_dir_tan <- file.path(curve_root, "Tangent")
if (isTRUE(save_mse_t)) {
  dir.create(curve_dir_lqd, recursive = TRUE, showWarnings = FALSE)
  dir.create(curve_dir_tan, recursive = TRUE, showWarnings = FALSE)
}

cat("\n===== baseline simulations (with optional MSE(t) curves) =====\n")
cat(sprintf("Start time: %s\n", as.character(Sys.time())))
cat(sprintf("reps=%d, seed_base=%d, t_design=%s\n", reps, seed_base, t_design))
cat(sprintf("grid: n=%s; N=%s; m=%s\n", paste(n_vec, collapse=","), paste(N_vec, collapse=","), paste(m_vec, collapse=",")))
cat(sprintf("Ku=%d, basis_type=%s, bs_degree=%d\n", Ku, basis_type, bs_degree))
cat(sprintf("smooth_deriv=%s, spar=%s\n", as.character(smooth_deriv), ifelse(is.null(spar), "NULL", as.character(spar))))
cat(sprintf("t_grid_len=%d, methodXi=%s\n", t_grid_len, methodXi))
cat(sprintf("decayR=%.3f (%s)\n", decayR,
            if (isTRUE(sim_supports_decayR)) "passed to simulate_transport_process" else "ignored (mod_script has no decayR arg)"))
cat(sprintf("cores=%d (SLURM_CPUS_PER_TASK=%s, SLURM_NTASKS=%s)\n",
            cores,
            ifelse(nchar(cpus_env)==0, "<unset>", cpus_env),
            ifelse(nchar(cores_env)==0, "<unset>", cores_env)))
cat(sprintf("outdir=%s\n", outdir))
cat(sprintf("mod_script=%s\n", mod_script))
cat(sprintf("base_script=%s\n", base_script))
cat(sprintf("save_mse_t=%s, plot_mse_t=%s, mse_rep_only=%s\n",
            as.character(save_mse_t), as.character(plot_mse_t),
            ifelse(is.na(mse_rep_only), "<all reps>", as.character(mse_rep_only))))
if (isTRUE(save_mse_t)) cat(sprintf("curve_root=%s\n", curve_root))
cat("============================================================\n\n")

# ---- helper: compute pointwise MSE(t) curve for a method ----
compute_mse_curve <- function(T_hat_list, T_true_list, t_grid, xgrid) {
  n <- length(T_hat_list)
  L <- length(t_grid)
  W1_sq <- matrix(NA_real_, nrow = n, ncol = L)
  for (i in seq_len(n)) {
    Th <- as.matrix(T_hat_list[[i]])
    Tt <- as.matrix(T_true_list[[i]])
    if (!all(dim(Th) == dim(Tt))) stop(sprintf("Subject %d: dim mismatch in T_hat vs T_true", i))
    # W1(t) for each row (time), then square
    w1_t <- trapz_rows(xgrid, abs(Th - Tt))
    W1_sq[i, ] <- w1_t^2
  }
  mse_t <- colMeans(W1_sq)
  se_t  <- if (n > 1) apply(W1_sq, 2, stats::sd) / sqrt(n) else rep(0, L)
  list(
    W1_sq = W1_sq,
    t = t_grid,
    mse = mse_t,
    se = se_t,
    lo = pmax(mse_t - 2 * se_t, 0),
    hi = mse_t + 2 * se_t
  )
}

write_curve_outputs <- function(curve, out_csv, out_pdf = NULL, main = "") {
  df <- data.frame(t = curve$t, mse = curve$mse, se = curve$se, lo = curve$lo, hi = curve$hi)
  utils::write.csv(df, out_csv, row.names = FALSE)

  if (!is.null(out_pdf)) {
    grDevices::pdf(out_pdf, width = 7, height = 4)
    graphics::plot(df$t, df$mse, type = "l", lwd = 2,
                   xlab = "t",
                   ylab = expression(paste("MSE(t) = mean ", d[W,1]^2, "(T_hat(t), T_true(t))")),
                   main = main)
    graphics::lines(df$t, df$lo, lty = 2)
    graphics::lines(df$t, df$hi, lty = 2)
    graphics::legend("topright", legend = c("MSE(t)", "± 2 SE"), lty = c(1, 2), bty = "n")
    graphics::grid()
    grDevices::dev.off()
  }
  invisible(df)
}

# ---- worker function for one replication ----
run_one_rep <- function(rep_id,
                        n_vec, N_vec, m_vec,
                        seed_base, t_design,
                        Ku, basis_type, bs_degree,
                        smooth_deriv, spar,
                        t_grid_out,
                        methodXi, decayR, sim_supports_decayR,beta_a,beta_b,
                        outdir, save_each_rep,
                        save_mse_t, plot_mse_t, mse_rep_only,
                        curve_dir_lqd, curve_dir_tan) {

  rows <- vector("list", length(n_vec) * length(N_vec) * length(m_vec))
  idx <- 0L

  do_save_curve_this_rep <- isTRUE(save_mse_t) && (is.na(mse_rep_only) || rep_id == mse_rep_only)

  for (n in n_vec) {
    for (N in N_vec) {
      for (m in m_vec) {

        seed <- seed_base + rep_id * 1000000L + n * 1000L + N * 10L + m

        out <- tryCatch({
          sim <- if (isTRUE(sim_supports_decayR)) {
            simulate_transport_process(
              n = n, N = N, m = m,
              t_design = t_design,
              seed = seed,
              decayR = decayR,
              beta_a=beta_a,
              beta_b=beta_b
            )
          } else {
            simulate_transport_process(
              n = n, N = N, m = m,
              t_design = t_design,
              seed = seed
            )
          }

          # Truth on common evaluation grid
          T_true <- lapply(seq_len(n), function(i) {
            true_transports_on_grid(
              Xi_row = sim$Xi[i, ],
              T0_row = sim$T0_true[i, ],
              t_grid = t_grid_out,
              xgrid  = sim$xgrid,
              K      = sim$K,
              g      = sim$g,
              project = TRUE
            )
          })

          # ---- LQD baseline
          fit_lqd <- fit_LQD_fpca(
            Q_list = sim$T_hat,
            t_list = sim$t_list,
            ugrid  = sim$xgrid,
            Ku = Ku,
            basis_type = basis_type,
            bs_degree = bs_degree,
            t_grid = t_grid_out,
            methodXi = methodXi,
            fpca_optns = list(methodMuCovEst = "smooth"),
            smooth_deriv = smooth_deriv,
            spar = spar,
            eps_deriv = 1e-8,
            lower = 0,
            upper = 1,
            project = TRUE,
            save_fpca_fits = FALSE,
            verbose = FALSE
          )

          imse_lqd_i <- imse_transports(fit_lqd$Q_hat, T_true, t_grid_out, sim$xgrid)
          imse_lqd <- mean(imse_lqd_i)

          # ---- Tangent/log-map baseline
          fit_tan <- fit_tangent_wass_fpca(
            T_list = sim$T_hat,
            t_list = sim$t_list,
            xgrid  = sim$xgrid,
            Ku = Ku,
            basis_type = basis_type,
            bs_degree = bs_degree,
            t_grid = t_grid_out,
            methodXi = methodXi,
            fpca_optns = list(methodMuCovEst = "smooth"),
            lower = 0,
            upper = 1,
            project = TRUE,
            save_fpca_fits = FALSE,
            verbose = FALSE
          )

          imse_tan_i <- imse_transports(fit_tan$T_hat, T_true, t_grid_out, sim$xgrid)
          imse_tan <- mean(imse_tan_i)

          # ---- Optional: pointwise MSE(t) curves (W1^2 averaged over subjects)
          if (do_save_curve_this_rep) {
            # Ensure dirs exist on worker
            dir.create(curve_dir_lqd, recursive = TRUE, showWarnings = FALSE)
            dir.create(curve_dir_tan, recursive = TRUE, showWarnings = FALSE)

            # LQD curve
            curve_lqd <- compute_mse_curve(fit_lqd$Q_hat, T_true, t_grid_out, sim$xgrid)
            out_csv_lqd <- file.path(curve_dir_lqd, sprintf("MSE_curve_rep%03d_n%d_N%d_m%d.csv", rep_id, n, N, m))
            out_pdf_lqd <- if (isTRUE(plot_mse_t)) {
              file.path(curve_dir_lqd, sprintf("MSE_curve_rep%03d_n%d_N%d_m%d.pdf", rep_id, n, N, m))
            } else NULL
            write_curve_outputs(
              curve_lqd,
              out_csv = out_csv_lqd,
              out_pdf = out_pdf_lqd,
              main = sprintf("LQD baseline: pointwise MSE(t)\n(rep=%d, n=%d, N=%d, m=%d)", rep_id, n, N, m)
            )

            # Tangent curve
            curve_tan <- compute_mse_curve(fit_tan$T_hat, T_true, t_grid_out, sim$xgrid)
            out_csv_tan <- file.path(curve_dir_tan, sprintf("MSE_curve_rep%03d_n%d_N%d_m%d.csv", rep_id, n, N, m))
            out_pdf_tan <- if (isTRUE(plot_mse_t)) {
              file.path(curve_dir_tan, sprintf("MSE_curve_rep%03d_n%d_N%d_m%d.pdf", rep_id, n, N, m))
            } else NULL
            write_curve_outputs(
              curve_tan,
              out_csv = out_csv_tan,
              out_pdf = out_pdf_tan,
              main = sprintf("Tangent baseline: pointwise MSE(t)\n(rep=%d, n=%d, N=%d, m=%d)", rep_id, n, N, m)
            )
          }

          list(ok = TRUE, imse_lqd = imse_lqd, imse_tan = imse_tan)

        }, error = function(e) {
          list(ok = FALSE, imse_lqd = NA_real_, imse_tan = NA_real_, err = conditionMessage(e))
        })

        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          rep = rep_id,
          n = n,
          N = N,
          m = m,
          imse_lqd = out$imse_lqd,
          imse_tangent = out$imse_tan,
          ok = out$ok,
          err = if (!isTRUE(out$ok)) out$err else "",
          stringsAsFactors = FALSE
        )

        rm(out)
        gc(verbose = FALSE)
      }
    }
  }

  df <- do.call(rbind, rows)

  if (isTRUE(save_each_rep)) {
    saveRDS(df, file = file.path(outdir, sprintf("rep_%03d.rds", rep_id)))
  }

  df
}

# ---- start cluster ----
cl <- makeCluster(cores)

# Export script paths
clusterExport(cl, varlist = c("mod_script", "base_script"), envir = environment())

# Load packages + source scripts on workers
invisible(clusterEvalQ(cl, {
  suppressMessages(library(fdapace))
  source(mod_script)
  source(base_script)
}))

clusterExport(
  cl,
  varlist = c("run_one_rep",
              "compute_mse_curve", "write_curve_outputs",
              "n_vec", "N_vec", "m_vec",
              "seed_base", "t_design",
              "Ku", "basis_type", "bs_degree",
              "smooth_deriv", "spar",
              "t_grid_out", "methodXi", "decayR", "sim_supports_decayR","beta_a","beta_b",
              "outdir", "save_each_rep",
              "save_mse_t", "plot_mse_t", "mse_rep_only",
              "curve_dir_lqd", "curve_dir_tan"),
  envir = environment()
)

rep_ids <- seq_len(reps)
cat(sprintf("Running %d replications in parallel...\n", reps))
ptm <- proc.time()
res_list <- parLapply(cl, rep_ids, function(r) {
  run_one_rep(r,
              n_vec, N_vec, m_vec,
              seed_base, t_design,
              Ku, basis_type, bs_degree,
              smooth_deriv, spar,
              t_grid_out,
              methodXi, decayR, sim_supports_decayR,beta_a,beta_b,
              outdir, save_each_rep,
              save_mse_t, plot_mse_t, mse_rep_only,
              curve_dir_lqd, curve_dir_tan)
})
run_time <- proc.time() - ptm

stopCluster(cl)

# ---- combine + save ----
raw_df <- do.call(rbind, res_list)
raw_path <- file.path(outdir, "raw_results.rds")
saveRDS(raw_df, raw_path)

# Summary by (n,N,m)
ok_df <- raw_df[raw_df$ok, , drop = FALSE]

summ_csv <- file.path(outdir, "summary.csv")

if (nrow(ok_df) == 0) {
  cat("\nWARNING: zero successful fits (ok_df has 0 rows).\n")
  cat("Writing raw_results.rds and an EMPTY summary.csv (check raw_results.rds$err).\n\n")

  err_tab <- sort(table(raw_df$err[raw_df$err != ""]), decreasing = TRUE)
  if (length(err_tab) > 0) {
    cat("Top error messages (up to 10):\n")
    print(utils::head(err_tab, 10))
    cat("\n")
  }

  summ_df <- data.frame()
  write.csv(summ_df, summ_csv, row.names = FALSE)

  sink(file.path(outdir, "sessionInfo.txt"))
  print(sessionInfo())
  sink()

  rdata_path <- file.path(outdir, "results_baselines.RData")
  save(
    raw_df, summ_df, ok_df, run_time,
    reps, seed_base, t_design,
    Ku, basis_type, bs_degree,
    smooth_deriv, spar,
    t_grid_len, t_grid_out,
    methodXi, decayR,
    n_vec, N_vec, m_vec,
    mod_script, base_script,
    save_mse_t, plot_mse_t, mse_rep_only,
    curve_root,
    sim_supports_decayR,
    file = rdata_path
  )

  cat("\n===== DONE =====\n")
  cat(sprintf("End time: %s\n", as.character(Sys.time())))
  cat(sprintf("Wall time (proc.time): %.2f seconds\n", as.numeric(run_time["elapsed"])))
  cat(sprintf("Raw results: %s\n", raw_path))
  cat(sprintf("Summary CSV: %s\n", summ_csv))
  if (isTRUE(save_mse_t)) {
    cat(sprintf("MSE(t) curves under: %s\n", curve_root))
    if (!is.na(mse_rep_only)) {
      cat(sprintf("(Only saved for rep %d; set mse_rep_only=NA to save for all reps.)\n", mse_rep_only))
    }
  }
  cat(sprintf("RData bundle: %s\n", rdata_path))
  cat("================\n")

  quit(status = 0)
}

split_key <- interaction(ok_df$n, ok_df$N, ok_df$m, drop = TRUE)
groups <- split(ok_df, split_key)

summ_list <- lapply(groups, function(d) {
  data.frame(
    n = d$n[1],
    N = d$N[1],
    m = d$m[1],
    reps_ok = nrow(d),
    mean_imse_lqd = mean(d$imse_lqd),
    se_imse_lqd = stats::sd(d$imse_lqd) / sqrt(nrow(d)),
    mean_imse_tangent = mean(d$imse_tangent),
    se_imse_tangent = stats::sd(d$imse_tangent) / sqrt(nrow(d)),
    stringsAsFactors = FALSE
  )
})

summ_df <- do.call(rbind, summ_list)
summ_df <- summ_df[order(summ_df$m, summ_df$n, summ_df$N), ]
write.csv(summ_df, summ_csv, row.names = FALSE)

sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()

# ---- save an .RData bundle ----
rdata_path <- file.path(outdir, "results_baselines.RData")
save(
  raw_df, summ_df, ok_df, run_time,
  reps, seed_base, t_design,
  Ku, basis_type, bs_degree,
  smooth_deriv, spar,
  t_grid_len, t_grid_out,
  methodXi,decayR,
  n_vec, N_vec, m_vec,
  mod_script, base_script,
  save_mse_t, plot_mse_t, mse_rep_only,
  curve_root,
  file = rdata_path
)

cat("\n===== DONE =====\n")
cat(sprintf("End time: %s\n", as.character(Sys.time())))
cat(sprintf("Wall time (proc.time): %.2f seconds\n", as.numeric(run_time["elapsed"])))
cat(sprintf("Raw results: %s\n", raw_path))
cat(sprintf("Summary CSV: %s\n", summ_csv))
if (isTRUE(save_mse_t)) {
  cat(sprintf("MSE(t) curves under: %s\n", curve_root))
  if (!is.na(mse_rep_only)) {
    cat(sprintf("(Only saved for rep %d; set mse_rep_only=NA to save for all reps.)\n", mse_rep_only))
  }
}
cat(sprintf("RData bundle: %s\n", rdata_path))
cat("================\n")
