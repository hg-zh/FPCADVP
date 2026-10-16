#!/usr/bin/env Rscript

## para_R2_mse_t_refined.R
## ------------------------------------------------------------
## Parallel Monte Carlo simulation driver for the PROPOSED method
## (transport-FPCA) implemented in transport_fpca_modular.R.
##
## This refined version fixes common issues:
##   - passes user-supplied latent/methodXi into fit_transport_fpca()
##   - supports the new simulation knob decayR
##   - computes IMSE on a FIXED evaluation grid t_grid_out (comparable to baselines)
##   - optionally saves pointwise MSE(t) curves per replication/setting
##   - avoids crashing when *all* fits fail (empty ok_df)
##   - robust boolean parsing (latent=F works)
##
## Outputs (in outdir):
##   - raw_results.rds
##   - summary.csv
##   - sessionInfo.txt
##   - (optional) MSE curves under outdir/mse_curves/Proposed/
##       MSE_curve_rep###_n##_N##_m##.csv
##       and (optional) matching PDFs if save_mse_plots=TRUE
##
## Example:
##   Rscript para_R2_mse_t_refined.R outdir=results/paper_comparison/proposed_latent_CE/setting2 \
##     reps=200 latent=TRUE methodXi=CE decayR=1 \
##     save_mse_t=TRUE save_mse_plots=FALSE t_grid_len=51 u_cap=0.95
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

parse_int_vec <- function(x) {
  x <- gsub("\\s+", "", as.character(x))
  if (nchar(x) == 0) return(integer(0))
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

parse_num_vec <- function(x) {
  s <- paste(as.character(x), collapse = ",")
  s <- gsub("\\s+", "", s)
  if (nchar(s) == 0) return(numeric(0))
  suppressWarnings(as.numeric(strsplit(s, ",", fixed = TRUE)[[1]]))
}


as_num_or_na <- function(x) {
  if (is.null(x)) return(NA_real_)
  s <- trimws(as.character(x)[1])
  if (tolower(s) %in% c("na","nan","inf","infinity","null","none","")) return(NA_real_)
  suppressWarnings(as.numeric(s))
}

as_int_or_na <- function(x) {
  if (is.null(x)) return(NA_integer_)
  s <- trimws(as.character(x)[1])
  if (tolower(s) %in% c("na","nan","inf","infinity","null","none","")) return(NA_integer_)
  suppressWarnings(as.integer(s))
}

as_method_or_null <- function(x) {
  if (is.null(x)) return(NULL)
  s <- trimws(as.character(x)[1])
  if (tolower(s) %in% c("", "auto", "null", "none")) return(NULL)
  s
}

call_with_supported_args <- function(fn, args) {
  fmls <- names(formals(fn))
  args2 <- args[names(args) %in% fmls]
  do.call(fn, args2)
}

# Default dependencies and outputs are relative to this distributed script.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this driver with Rscript.")
script_dir <- dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg), fixed = TRUE), mustWork = TRUE))

kv <- parse_kv_args(commandArgs(trailingOnly = TRUE))

# ---- config (defaults) ----
reps       <- if (!is.null(kv$reps)) as.integer(kv$reps) else 200
seed_base  <- if (!is.null(kv$seed_base)) as.integer(kv$seed_base) else 2000
outdir     <- kv$outdir %||% file.path(script_dir, "results", "paper_comparison", "proposed_latent_CE", "setting2")
t_design   <- kv$t_design %||% "random"   # random|fixed

latent     <- if (!is.null(kv$latent)) as_bool(kv$latent, default = TRUE) else TRUE
methodXi   <- as_method_or_null(kv$methodXi)  # NULL lets fit_transport_fpca choose

# NEW simulation knob
decayR     <- if (!is.null(kv$decayR)) as.numeric(kv$decayR) else 1
if (!is.finite(decayR) || decayR <= 0) decayR <- 1

# Proposed-method stability knobs
eps_u <- if (!is.null(kv$eps_u)) as_num_or_na(kv$eps_u) else 0.05
u_cap <- if (!is.null(kv$u_cap)) as_num_or_na(kv$u_cap) else NA_real_
Z_cap <- if (!is.null(kv$Z_cap)) as_num_or_na(kv$Z_cap) else 8

beta_a <- if (!is.null(kv$beta_a) && nchar(trimws(kv$beta_a)) > 0) parse_num_vec(kv$beta_a) else c(3,4)
beta_b <- if (!is.null(kv$beta_b) && nchar(trimws(kv$beta_b)) > 0) parse_num_vec(kv$beta_b) else c(1,2)

for (beta_name in c("beta_a", "beta_b")) {
  beta_range <- get(beta_name)
  if (length(beta_range) != 2 || any(!is.finite(beta_range)) ||
      any(beta_range <= 0) || beta_range[1] >= beta_range[2]) {
    stop(beta_name, " must be a finite, positive, increasing pair")
  }
}

save_each_rep <- if (!is.null(kv$save_each_rep)) as_bool(kv$save_each_rep, TRUE) else TRUE

# ---- MSE(t) curve options ----
save_mse_t      <- if (!is.null(kv$save_mse_t)) as_bool(kv$save_mse_t, FALSE) else FALSE
save_mse_plots  <- if (!is.null(kv$save_mse_plots)) as_bool(kv$save_mse_plots, FALSE) else FALSE
mse_rep_only    <- as_int_or_na(kv$mse_rep_only)
t_grid_len      <- if (!is.null(kv$t_grid_len)) as.integer(kv$t_grid_len) else 51
if (!is.finite(t_grid_len) || t_grid_len < 2) t_grid_len <- 51
t_grid_out <- seq(0, 1, length.out = t_grid_len)

# Parameter grid (defaults match your current baseline grid)
n_vec <- if (!is.null(kv$n_vec)) parse_int_vec(kv$n_vec) else c(20, 50, 100, 200, 400)
N_vec <- if (!is.null(kv$N_vec)) parse_int_vec(kv$N_vec) else c(3, 5, 10, 20)
m_vec <- if (!is.null(kv$m_vec)) parse_int_vec(kv$m_vec) else c(200)

if (length(n_vec) < 1 || any(!is.finite(n_vec))) stop("n_vec is empty or invalid")
if (length(N_vec) < 1 || any(!is.finite(N_vec))) stop("N_vec is empty or invalid")
if (length(m_vec) < 1 || any(!is.finite(m_vec))) stop("m_vec is empty or invalid")

# ---- cores: prefer CPUS_PER_TASK for shared-memory jobs ----
cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
if (!is.finite(cores) || cores < 1) cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_NTASKS")))
if ((!is.finite(cores) || cores < 1) && !is.null(kv$cores)) cores <- as.integer(kv$cores)
if (!is.finite(cores) || cores < 1) cores <- 1

cores_env <- Sys.getenv("SLURM_NTASKS")
cpus_env  <- Sys.getenv("SLURM_CPUS_PER_TASK")

# ---- locate modular script ----
mod_script <- kv$mod_script %||% file.path(script_dir, "transport_fpca_modular.R")
if (!file.exists(mod_script)) {
  stop(sprintf("Cannot find '%s'. Put transport_fpca_modular.R in the same directory or pass mod_script=PATH.", mod_script))
}
mod_script <- normalizePath(mod_script)

# ---- preflight check: ensure mod_script supports decayR and key args ----
tmp_env <- new.env(parent = globalenv())
sys.source(mod_script, envir = tmp_env)
if (!exists("simulate_transport_process", envir = tmp_env, inherits = FALSE)) {
  stop("mod_script does not define simulate_transport_process().")
}
if (!("decayR" %in% names(formals(get("simulate_transport_process", envir = tmp_env))))) {
  stop("Your transport_fpca_modular.R does not support the argument decayR.
",
       "Use the updated modular file (the one used by baselines) or switch mod_script to the correct path.")
}
if (!exists("fit_transport_fpca", envir = tmp_env, inherits = FALSE)) {
  stop("mod_script does not define fit_transport_fpca().")
}
need_args <- c("latent", "methodXi")
have_args <- names(formals(get("fit_transport_fpca", envir = tmp_env)))
if (!all(need_args %in% have_args)) {
  stop("fit_transport_fpca() in mod_script is missing expected args: ",
       paste(setdiff(need_args, have_args), collapse = ", "))
}

rm(tmp_env)


# ---- output dirs ----
if (file.exists(file.path(outdir, "raw_results.rds"))) stop("Refusing to overwrite existing results; choose a new outdir: ", outdir)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

curve_root <- kv$curve_root %||% file.path(outdir, "mse_curves")
curve_dir_prop <- file.path(curve_root, "Proposed")
if (isTRUE(save_mse_t)) {
  dir.create(curve_dir_prop, recursive = TRUE, showWarnings = FALSE)
}

cat("\n===== transport-FPCA parallel simulation =====\n")
cat(sprintf("Start time: %s\n", as.character(Sys.time())))
cat(sprintf("reps=%d, seed_base=%d, t_design=%s, latent=%s, methodXi=%s\n",
            reps, seed_base, t_design,
            as.character(latent),
            ifelse(is.null(methodXi), "<AUTO>", methodXi)))
cat(sprintf("grid: n=%s; N=%s; m=%s\n", paste(n_vec, collapse=","), paste(N_vec, collapse=","), paste(m_vec, collapse=",")))
cat(sprintf("cores=%d (SLURM_NTASKS=%s, SLURM_CPUS_PER_TASK=%s)\n",
            cores,
            ifelse(nchar(cores_env)==0, "<unset>", cores_env),
            ifelse(nchar(cpus_env)==0, "<unset>", cpus_env)))
cat(sprintf("outdir=%s\n", outdir))
cat(sprintf("mod_script=%s\n", mod_script))
cat(sprintf("decayR=%.3f\n", decayR))
cat(sprintf("save_mse_t=%s, save_mse_plots=%s, t_grid_len=%d, mse_rep_only=%s\n",
            as.character(save_mse_t),
            as.character(save_mse_plots),
            t_grid_len,
            ifelse(is.na(mse_rep_only), "<all reps>", as.character(mse_rep_only))))
cat(sprintf("eps_u=%.4f, u_cap=%s, Z_cap=%s\n",
            as.numeric(eps_u),
            ifelse(is.finite(u_cap), sprintf("%.4f", u_cap), "<AUTO>"),
            ifelse(is.finite(Z_cap), sprintf("%.4f", Z_cap), "<NULL>")))
cat("===========================================\n\n")

# ---- helpers: interpolate a transport matrix in time ----
interp_Tmat_time <- function(Th, t_from, t_to) {
  Th <- as.matrix(Th)
  t_from <- as.numeric(t_from)
  t_to <- as.numeric(t_to)
  G <- ncol(Th)
  out <- matrix(NA_real_, nrow = length(t_to), ncol = G)
  for (gg in seq_len(G)) {
    out[, gg] <- approx(x = t_from, y = Th[, gg], xout = t_to, rule = 2)$y
  }
  out
}

compute_mse_curve <- function(T_hat_list, T_true_list, t_grid, xgrid) {
  n <- length(T_hat_list)
  L <- length(t_grid)
  W1_sq <- matrix(NA_real_, nrow = n, ncol = L)
  for (i in seq_len(n)) {
    Th <- as.matrix(T_hat_list[[i]])
    Tt <- as.matrix(T_true_list[[i]])
    if (!all(dim(Th) == dim(Tt))) stop(sprintf("Subject %d: dim mismatch in T_hat vs T_true", i))
    w1_t <- trapz_rows(xgrid, abs(Th - Tt))
    W1_sq[i, ] <- w1_t^2
  }
  mse_t <- colMeans(W1_sq)
  se_t  <- if (n > 1) apply(W1_sq, 2, stats::sd) / sqrt(n) else rep(0, L)
  list(
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
                        latent, methodXi,
                        eps_u, u_cap, Z_cap,
                        t_grid_out,
                        decayR,beta_a,beta_b,
                        outdir, save_each_rep,
                        save_mse_t, save_mse_plots, mse_rep_only,
                        curve_dir_prop) {

  rows <- vector("list", length(n_vec) * length(N_vec) * length(m_vec))
  idx <- 0L

  do_save_curve_this_rep <- isTRUE(save_mse_t) && (is.na(mse_rep_only) || rep_id == mse_rep_only)

  for (n in n_vec) {
    for (N in N_vec) {
      for (m in m_vec) {

        seed <- seed_base + rep_id * 1000000L + n * 1000L + N * 10L + m

        out <- tryCatch({
          sim <- simulate_transport_process(
            n = n, N = N, m = m,
            t_design = t_design,
            seed = seed,
            decayR = decayR,
            beta_a = beta_a,
            beta_b = beta_b
          )

          # Fit proposed method
          fit_args <- list(
            T_list = sim$T_hat,
            t_list = sim$t_list,
            xgrid  = sim$xgrid,
            kappa  = sim$kappa_true,
            latent = latent,
            methodXi = methodXi,
            fpca_optns = list(methodMuCovEst = "smooth"),
            eps_u = eps_u,
            project = TRUE,
            verbose = FALSE
          )
          if (is.finite(u_cap)) fit_args$u_cap <- u_cap
          if (is.finite(Z_cap)) fit_args$Z_cap <- Z_cap

          fit <- call_with_supported_args(fit_transport_fpca, fit_args)

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

          # Put T_hat onto t_grid_out if needed
          T_hat <- fit$T_hat
          if (length(fit$workGrid) != length(t_grid_out) || max(abs(fit$workGrid - t_grid_out)) > 1e-12) {
            T_hat <- lapply(T_hat, function(Th) interp_Tmat_time(Th, fit$workGrid, t_grid_out))
          }

          # Pointwise W1(t) per subject, then integrate over t
          # (this matches the baseline IMSE definition)
          imse_i <- numeric(n)
          imse_sq_i <- numeric(n)
          for (i in seq_len(n)) {
            Th <- as.matrix(T_hat[[i]])
            Tt <- as.matrix(T_true[[i]])
            if (!all(dim(Th) == dim(Tt))) stop("dim mismatch after interpolation")
            w1_t <- trapz_rows(sim$xgrid, abs(Th - Tt))
            imse_i[i] <- trapz_vec(t_grid_out, w1_t)
            imse_sq_i[i] <- trapz_vec(t_grid_out, w1_t^2)
          }

          # Save MSE(t) curve for this replication/setting (optional)
          if (do_save_curve_this_rep) {
            dir.create(curve_dir_prop, recursive = TRUE, showWarnings = FALSE)
            curve <- compute_mse_curve(T_hat, T_true, t_grid_out, sim$xgrid)
            out_csv <- file.path(curve_dir_prop, sprintf("MSE_curve_rep%03d_n%d_N%d_m%d.csv", rep_id, n, N, m))
            out_pdf <- if (isTRUE(save_mse_plots)) {
              file.path(curve_dir_prop, sprintf("MSE_curve_rep%03d_n%d_N%d_m%d.pdf", rep_id, n, N, m))
            } else NULL
            write_curve_outputs(
              curve,
              out_csv = out_csv,
              out_pdf = out_pdf,
              main = sprintf("Proposed: pointwise MSE(t)\n(rep=%d, n=%d, N=%d, m=%d)", rep_id, n, N, m)
            )
          }

          list(ok = TRUE,
               imse = mean(imse_i),
               imse_sq = mean(imse_sq_i))

        }, error = function(e) {
          list(ok = FALSE, imse = NA_real_, imse_sq = NA_real_, err = conditionMessage(e))
        })

        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          rep = rep_id,
          n = n,
          N = N,
          m = m,
          imse = out$imse,
          imse_sq = out$imse_sq,
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

clusterExport(cl, varlist = c("mod_script"), envir = environment())

invisible(clusterEvalQ(cl, {
  suppressMessages(library(fdapace))
  source(mod_script)
}))

clusterExport(
  cl,
  varlist = c("run_one_rep",
              "interp_Tmat_time", "compute_mse_curve", "write_curve_outputs",
              "call_with_supported_args",
              "n_vec", "N_vec", "m_vec",
              "seed_base", "t_design",
              "latent", "methodXi",
              "eps_u", "u_cap", "Z_cap",
              "t_grid_out",
              "decayR","beta_a","beta_b",
              "outdir", "save_each_rep",
              "save_mse_t", "save_mse_plots", "mse_rep_only",
              "curve_dir_prop"),
  envir = environment()
)

rep_ids <- seq_len(reps)
cat(sprintf("Running %d replications in parallel...\n", reps))
ptm <- proc.time()
res_list <- parLapply(cl, rep_ids, function(r) {
  run_one_rep(r,
              n_vec, N_vec, m_vec,
              seed_base, t_design,
              latent, methodXi,
              eps_u, u_cap, Z_cap,
              t_grid_out,
              decayR,beta_a,beta_b,
              outdir, save_each_rep,
              save_mse_t, save_mse_plots, mse_rep_only,
              curve_dir_prop)
})
run_time <- proc.time() - ptm

stopCluster(cl)

# ---- combine + save ----
raw_df <- do.call(rbind, res_list)
raw_path <- file.path(outdir, "raw_results.rds")
saveRDS(raw_df, raw_path)

# Summary: mean + Monte Carlo SE by (n,N,m)
ok_df <- raw_df[raw_df$ok, , drop = FALSE]

summ_csv <- file.path(outdir, "summary.csv")

if (nrow(ok_df) == 0) {
  cat("\nWARNING: zero successful fits (ok_df has 0 rows).\n")
  cat("Writing raw_results.rds and an EMPTY summary.csv (check raw_results.rds$err).\n")
  write.csv(data.frame(), summ_csv, row.names = FALSE)

  # Print a quick summary of the failure reasons
  err_vec <- raw_df$err
  err_vec <- err_vec[nzchar(err_vec)]
  if (length(err_vec) > 0) {
    tab <- sort(table(err_vec), decreasing = TRUE)
    cat("\nTop error messages (up to 10):\n")
    print(head(tab, 10))
  }
} else {
  split_key <- interaction(ok_df$n, ok_df$N, ok_df$m, drop = TRUE)
  groups <- split(ok_df, split_key)

  summ_list <- lapply(groups, function(d) {
    reps_ok <- nrow(d)
    se_imse <- if (reps_ok > 1) stats::sd(d$imse) / sqrt(reps_ok) else 0
    se_imse_sq <- if (reps_ok > 1) stats::sd(d$imse_sq) / sqrt(reps_ok) else 0
    data.frame(
      n = d$n[1],
      N = d$N[1],
      m = d$m[1],
      reps_ok = reps_ok,
      mean_imse = mean(d$imse),
      se_imse = se_imse,
      mean_imse_sq = mean(d$imse_sq),
      se_imse_sq = se_imse_sq,
      stringsAsFactors = FALSE
    )
  })

  summ_df <- do.call(rbind, summ_list)
  summ_df <- summ_df[order(summ_df$m, summ_df$n, summ_df$N), ]
  write.csv(summ_df, summ_csv, row.names = FALSE)
}

# Save session info
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()

# ---- save everything to .RData ----
rdata_path <- file.path(outdir, "results_prop.RData")
save(
  raw_df, ok_df,
  reps, seed_base, t_design,
  latent, methodXi,
  decayR,
  eps_u, u_cap, Z_cap,
  n_vec, N_vec, m_vec,
  beta_a, beta_b,
  t_grid_len, t_grid_out,
  cores, run_time, mod_script,
  save_mse_t, save_mse_plots, mse_rep_only,
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
cat("===============\n")
