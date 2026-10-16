# transport_fpca_modular_robust.R
# ------------------------------------------------------------
# Robustified version of transport_fpca_modular.R.
#
# Key change: protect the latent-Gaussian (sparse) pipeline from
# blow-ups caused by g^{-1}(u) when u is very close to +/-1.
#
# In simulations with estimated transports, noise can yield |U_hat|>1.
# If you clip to 1-eps with eps extremely small, then g^{-1} produces
# enormous values that dominate FPCA, and this gets worse as n grows.
#
# This version adds:
#   - u_cap: clamp U_hat to [-u_cap, u_cap] before applying g^{-1}
#   - Z_cap: optional winsorization cap on Z=g^{-1}(U)
#   - diagnostics returned in fit$diagnostics
# ------------------------------------------------------------

# =========================
# 0) Small helpers
# =========================

clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

trapz_vec <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) != length(y)) stop("trapz_vec: length mismatch")
  if (length(x) < 2) return(0)
  dx <- diff(x)
  sum((y[-length(y)] + y[-1]) * dx / 2)
}

trapz_rows <- function(x, Y) {
  x <- as.numeric(x)
  Y <- as.matrix(Y)
  G <- length(x)
  if (ncol(Y) != G) stop("trapz_rows: ncol(Y) must equal length(x)")
  if (G < 2) return(rep(0, nrow(Y)))
  dx <- diff(x)
  mids <- (Y[, -G, drop = FALSE] + Y[, -1, drop = FALSE]) / 2
  mids <- sweep(mids, 2, dx, "*")
  rowSums(mids)
}

project_transport <- function(T, lower = 0, upper = 1) {
  T <- as.numeric(T)
  if (length(T) < 2) return(T)
  T <- clamp(T, lower, upper)
  T[1] <- lower
  T[length(T)] <- upper
  T <- cummax(T)
  T[length(T)] <- upper
  T[1] <- lower
  T
}

inverse_transport <- function(T, xgrid) {
  T <- as.numeric(T)
  xgrid <- as.numeric(xgrid)
  if (length(T) != length(xgrid)) stop("inverse_transport: length mismatch")
  Tm <- project_transport(T, lower = min(xgrid), upper = max(xgrid))
  idx <- !duplicated(Tm)
  x_u <- Tm[idx]
  y_u <- xgrid[idx]
  if (length(x_u) < 2) return(xgrid)
  approx(x = x_u, y = y_u, xout = xgrid, rule = 2, ties = "ordered")$y
}

transport_sign_rows <- function(Tmat, xgrid, tol = 1e-12) {
  xgrid <- as.numeric(xgrid)
  Tmat <- as.matrix(Tmat)
  if (ncol(Tmat) != length(xgrid)) stop("transport_sign_rows: grid mismatch")
  disp <- sweep(Tmat, 2, xgrid, "-")
  ints <- trapz_rows(xgrid, disp)
  s <- sign(ints)
  s[abs(ints) < tol] <- 1
  s
}

transport_norm1_rows <- function(Tmat, xgrid) {
  xgrid <- as.numeric(xgrid)
  Tmat <- as.matrix(Tmat)
  disp <- abs(sweep(Tmat, 2, xgrid, "-"))
  trapz_rows(xgrid, disp)
}

w1_distance_vec <- function(T1, T2, xgrid) {
  xgrid <- as.numeric(xgrid)
  trapz_vec(xgrid, abs(as.numeric(T1) - as.numeric(T2)))
}

rescale_transport_to_norm <- function(T, xgrid, target_norm, project = TRUE, tol = 1e-12) {
  xgrid <- as.numeric(xgrid)
  T <- as.numeric(T)
  cur <- trapz_vec(xgrid, abs(T - xgrid))
  if (!is.finite(cur) || cur < tol) {
    out <- xgrid
    return(if (project) project_transport(out) else out)
  }
  out <- xgrid + (target_norm / cur) * (T - xgrid)
  if (project) out <- project_transport(out)
  out
}

scalar_mult_transport <- function(alpha, T0, xgrid, project = TRUE) {
  xgrid <- as.numeric(xgrid)
  T0 <- as.numeric(T0)
  if (length(T0) != length(xgrid)) stop("scalar_mult_transport: grid mismatch")
  if (!is.finite(alpha)) stop("scalar_mult_transport: alpha must be finite")

  if (alpha >= 0) {
    out <- xgrid + alpha * (T0 - xgrid)
  } else {
    T0_inv <- inverse_transport(T0, xgrid)
    out <- xgrid + (-alpha) * (T0_inv - xgrid)
  }
  if (project) out <- project_transport(out)
  out
}

# =========================
# 1) Baseline estimation per subject
# =========================

estimate_subject_baselines <- function(T_i, xgrid, kappa, project = TRUE, tol_sign = 1e-12) {
  T_i <- as.matrix(T_i)
  signs <- transport_sign_rows(T_i, xgrid, tol = tol_sign)
  I_pos <- which(signs > 0)
  I_neg <- which(signs < 0)

  T_plus_raw  <- if (length(I_pos) > 0) colMeans(T_i[I_pos, , drop = FALSE]) else rep(NA_real_, length(xgrid))
  T_minus_raw <- if (length(I_neg) > 0) colMeans(T_i[I_neg, , drop = FALSE]) else rep(NA_real_, length(xgrid))

  T_plus  <- if (!anyNA(T_plus_raw))  rescale_transport_to_norm(T_plus_raw,  xgrid, kappa, project = project) else rep(NA_real_, length(xgrid))
  T_minus <- if (!anyNA(T_minus_raw)) rescale_transport_to_norm(T_minus_raw, xgrid, kappa, project = project) else rep(NA_real_, length(xgrid))

  if (anyNA(T_plus) && !anyNA(T_minus))  T_plus  <- inverse_transport(T_minus, xgrid)
  if (anyNA(T_minus) && !anyNA(T_plus))  T_minus <- inverse_transport(T_plus, xgrid)

  list(
    T_plus  = T_plus,
    T_minus = T_minus,
    signs   = signs,
    I_pos   = I_pos,
    I_neg   = I_neg
  )
}

# =========================
# 2) Main method: transport FPCA (robust)
# =========================

fit_transport_fpca <- function(T_list,
                               t_list,
                               xgrid,
                               kappa = NULL,
                               latent = TRUE,
                               g = function(z) 2 * atan(z) / pi,
                               g_inv = function(u) tan(pi * u / 2),
                               fpca_optns = NULL,
                               methodXi = NULL,
                               # Two ways to keep U away from +/-1:
                               #   - set u_cap directly (recommended)
                               #   - or use eps_u, where u_cap = 1 - eps_u
                               u_cap = NULL,
                               eps_u = 0.05,
                               # Optional cap on Z = g^{-1}(U)
                               Z_cap = 8,
                               project = TRUE,
                               tol_sign = 1e-12,
                               verbose = FALSE) {

  if (!requireNamespace("fdapace", quietly = TRUE)) {
    stop("Package 'fdapace' is required. Install it via install.packages('fdapace').")
  }

  if (!is.list(T_list) || !is.list(t_list)) stop("T_list and t_list must be lists.")
  n <- length(T_list)
  if (length(t_list) != n) stop("Length mismatch: T_list and t_list must have same length.")
  xgrid <- as.numeric(xgrid)
  G <- length(xgrid)

  # Determine u_cap
  if (!is.null(u_cap)) {
    u_cap <- as.numeric(u_cap)
    if (!is.finite(u_cap) || u_cap <= 0 || u_cap >= 1) stop("u_cap must be in (0,1).")
  } else {
    eps_u <- as.numeric(eps_u)
    if (!is.finite(eps_u) || eps_u <= 0 || eps_u >= 1) stop("eps_u must be in (0,1).")
    u_cap <- 1 - eps_u
  }

  # Reorder each subject by time
  for (i in seq_len(n)) {
    ti <- as.numeric(t_list[[i]])
    Ti <- as.matrix(T_list[[i]])
    if (ncol(Ti) != G) stop(sprintf("Subject %d: ncol(T_i) must equal length(xgrid).", i))
    if (length(ti) != nrow(Ti)) stop(sprintf("Subject %d: length(t_i) must equal nrow(T_i).", i))
    o <- order(ti)
    t_list[[i]] <- ti[o]
    T_list[[i]] <- Ti[o, , drop = FALSE]
  }

  # Choose kappa if not supplied
  if (is.null(kappa)) {
    all_norms <- unlist(lapply(T_list, function(Ti) transport_norm1_rows(Ti, xgrid)))
    kappa <- max(all_norms, na.rm = TRUE)
    if (!is.finite(kappa) || kappa <= 0) stop("Could not infer a valid kappa from the data.")
    if (verbose) message(sprintf("kappa not supplied; using max observed norm: kappa = %.4f", kappa))
  }

  # 1) Estimate baselines
  T0_plus  <- matrix(NA_real_, nrow = n, ncol = G)
  T0_minus <- matrix(NA_real_, nrow = n, ncol = G)
  signs_list <- vector("list", n)

  for (i in seq_len(n)) {
    base_i <- estimate_subject_baselines(T_list[[i]], xgrid, kappa, project = project, tol_sign = tol_sign)
    T0_plus[i, ]  <- base_i$T_plus
    T0_minus[i, ] <- base_i$T_minus
    signs_list[[i]] <- base_i$signs
  }

  # 2) Estimate U_hat at observed times, KEEP both raw and clamped
  U_raw_list <- vector("list", n)
  U_hat_list <- vector("list", n)

  for (i in seq_len(n)) {
    norms <- transport_norm1_rows(T_list[[i]], xgrid)
    u_raw <- as.numeric(signs_list[[i]]) * norms / kappa
    U_raw_list[[i]] <- u_raw
    U_hat_list[[i]] <- clamp(u_raw, -u_cap, u_cap)
  }

  # Diagnostics: how often do we exceed [-1,1] before clamping?
  u_raw_all <- unlist(U_raw_list)
  diagnostics <- list(
    n_subjects = n,
    n_obs_total = sum(lengths(t_list)),
    kappa = kappa,
    u_cap = u_cap,
    max_abs_U_raw = max(abs(u_raw_all)),
    frac_abs_U_raw_gt1 = mean(abs(u_raw_all) > 1),
    max_abs_U_used = max(abs(unlist(U_hat_list)))
  )

  # 3) FPCA on Z=g^{-1}(U) (latent) or directly on U (dense-style)
  if (latent) {
    Ly <- lapply(U_hat_list, g_inv)
    if (!is.null(Z_cap)) {
      Z_cap <- as.numeric(Z_cap)
      if (!is.finite(Z_cap) || Z_cap <= 0) stop("Z_cap must be a positive number or NULL")
      Ly <- lapply(Ly, function(z) clamp(z, -Z_cap, Z_cap))
      diagnostics$Z_cap <- Z_cap
    } else {
      diagnostics$Z_cap <- NA_real_
    }
    diagnostics$max_abs_Z_used <- max(abs(unlist(Ly)))
  } else {
    Ly <- U_hat_list
    diagnostics$Z_cap <- NA_real_
    diagnostics$max_abs_Z_used <- NA_real_
  }
  Lt <- t_list

  # 4) Fit FPCA
  if (is.null(methodXi)) methodXi <- if (latent) "CE" else "IN"
  if (is.null(fpca_optns)) fpca_optns <- list()
  fpca_optns$methodXi <- methodXi
  if (is.null(fpca_optns$methodMuCovEst)) fpca_optns$methodMuCovEst <- "smooth"

  fpca_fit <- fdapace::FPCA(Ly, Lt, optns = fpca_optns)

  # 5) Reconstruct scalar trajectories
  workGrid <- as.numeric(fpca_fit$workGrid)
  L <- length(workGrid)
  mu <- as.numeric(fpca_fit$mu)
  if (length(mu) != L) {
    mu <- as.numeric(mu)
    if (length(mu) != L) stop("fdapace::FPCA returned mu with unexpected length.")
  }

  traj_hat <- fpca_fit$xiEst %*% t(fpca_fit$phi)
  traj_hat <- sweep(traj_hat, 2, mu, "+")

  if (latent) {
    Z_hat_grid <- traj_hat
    if (!is.null(Z_cap)) Z_hat_grid <- clamp(Z_hat_grid, -Z_cap, Z_cap)
    U_hat_grid <- g(Z_hat_grid)
  } else {
    Z_hat_grid <- NULL
    U_hat_grid <- clamp(traj_hat, -u_cap, u_cap)
  }

  # 6) Reconstruct transports
  T_hat_list <- vector("list", n)
  for (i in seq_len(n)) {
    Ti_hat <- matrix(NA_real_, nrow = L, ncol = G)
    for (ell in seq_len(L)) {
      u_val <- U_hat_grid[i, ell]
      if (u_val >= 0) {
        Ti_hat[ell, ] <- xgrid + u_val * (T0_plus[i, ] - xgrid)
      } else {
        Ti_hat[ell, ] <- xgrid + (-u_val) * (T0_minus[i, ] - xgrid)
      }
      if (project) Ti_hat[ell, ] <- project_transport(Ti_hat[ell, ])
    }
    T_hat_list[[i]] <- Ti_hat
  }

  structure(list(
    xgrid = xgrid,
    kappa = kappa,
    latent = latent,
    g = g,
    g_inv = g_inv,
    project = project,
    u_cap = u_cap,
    eps_u = 1 - u_cap,
    Z_cap = if (is.null(Z_cap)) NA_real_ else Z_cap,
    diagnostics = diagnostics,
    t_list = t_list,
    T_obs = T_list,
    signs_obs = signs_list,
    T0_plus = T0_plus,
    T0_minus = T0_minus,
    U_hat_raw = U_raw_list,
    U_hat_obs = U_hat_list,
    fpca_fit = fpca_fit,
    workGrid = workGrid,
    Z_hat = Z_hat_grid,
    U_hat = U_hat_grid,
    T_hat = T_hat_list,
    call = match.call()
  ), class = "transport_fpca_fit")
}

# =========================
# 3) Simulation utilities (unchanged)
# =========================

create_cosine_basis <- function(t, K) {
  t <- as.numeric(t)
  if (K < 1) stop("K must be >= 1")
  Phi <- matrix(0, nrow = length(t), ncol = K)
  Phi[, 1] <- 1
  if (K >= 2) {
    Phi[, 2:K] <- cos(2 * pi * outer(t, 1:(K - 1)))
  }
  Phi
}

sample_from_transport <- function(T, xgrid, m) {
  u <- runif(m, 0, 1)
  approx(x = xgrid, y = as.numeric(T), xout = u, rule = 2)$y
}

estimate_transport_from_samples <- function(samples, probs, enforce_bounds = TRUE) {
  q <- as.numeric(stats::quantile(samples, probs = probs, type = 8, names = FALSE, na.rm = TRUE))
  if (enforce_bounds) {
    q[1] <- 0
    q[length(q)] <- 1
    q <- project_transport(q)
  }
  q
}

estimate_transports_from_sample_matrix <- function(sample_mat, probs, enforce_bounds = TRUE) {
  sample_mat <- as.matrix(sample_mat)
  Ni <- nrow(sample_mat)
  G <- length(probs)
  out <- matrix(NA_real_, nrow = Ni, ncol = G)
  for (j in seq_len(Ni)) {
    out[j, ] <- estimate_transport_from_samples(sample_mat[j, ], probs, enforce_bounds = enforce_bounds)
  }
  out
}

simulate_transport_process <- function(n,
                                       N,
                                       m,
                                       decayR = 1,
                                       xgrid = seq(0, 1, length.out = 128),
                                       K = 50,
                                       t_design = c("random", "fixed"),
                                       g = function(z) 2 * atan(z) / pi,
                                       beta_a = c(3, 4),
                                       beta_b = c(1, 2),
                                       seed = NULL,
                                       project = TRUE) {

  t_design <- match.arg(t_design)
  decayR <- as.numeric(decayR)
  if (!is.finite(decayR) || decayR <= 0) stop("simulate_transport_process: decayR must be a positive number")
  if (!is.null(seed)) set.seed(seed)

  xgrid <- as.numeric(xgrid)
  G <- length(xgrid)

  t_list <- vector("list", n)
  if (t_design == "random") {
    for (i in seq_len(n)) t_list[[i]] <- sort(runif(N, 0, 1))
  } else {
    tg <- seq(0, 1, length.out = N + 2)[2:(N + 1)]
    for (i in seq_len(n)) t_list[[i]] <- tg
  }

  Xi <- sapply(1:K, function(k) rnorm(n, mean = 0, sd = 1 / (k^decayR)))
  if (!is.matrix(Xi)) Xi <- matrix(Xi, nrow = n, ncol = K)

  a <- runif(n, beta_a[1], beta_a[2])
  b <- runif(n, beta_b[1], beta_b[2])
  T0_raw <- matrix(NA_real_, nrow = n, ncol = G)
  for (i in seq_len(n)) T0_raw[i, ] <- stats::qbeta(xgrid, a[i], b[i])

  norms0 <- transport_norm1_rows(T0_raw, xgrid)
  kappa_true <- mean(norms0)
  T0 <- matrix(NA_real_, nrow = n, ncol = G)
  for (i in seq_len(n)) T0[i, ] <- rescale_transport_to_norm(T0_raw[i, ], xgrid, kappa_true, project = project)

  Z_obs <- vector("list", n)
  U_obs <- vector("list", n)
  T_true <- vector("list", n)
  X_samp <- vector("list", n)
  T_hat <- vector("list", n)

  for (i in seq_len(n)) {
    Phi_i <- create_cosine_basis(t_list[[i]], K)
    Zi <- as.numeric(Phi_i %*% Xi[i, ])
    Ui <- g(Zi)

    Z_obs[[i]] <- Zi
    U_obs[[i]] <- Ui

    Ti_true <- matrix(NA_real_, nrow = N, ncol = G)
    Xi_mat <- matrix(NA_real_, nrow = N, ncol = m)

    for (j in seq_len(N)) {
      Ti_true[j, ] <- scalar_mult_transport(Ui[j], T0[i, ], xgrid, project = project)
      Xi_mat[j, ] <- sample_from_transport(Ti_true[j, ], xgrid, m)
    }

    T_true[[i]] <- Ti_true
    X_samp[[i]] <- Xi_mat
    T_hat[[i]] <- estimate_transports_from_sample_matrix(Xi_mat, xgrid, enforce_bounds = TRUE)
  }

  list(
    xgrid = xgrid,
    t_list = t_list,
    K = K,
    Xi = Xi,
    g = g,
    kappa_true = kappa_true,
    T0_true = T0,
    Z_obs = Z_obs,
    U_obs = U_obs,
    T_true = T_true,
    X_samp = X_samp,
    T_hat = T_hat
  )
}

true_transports_on_grid <- function(Xi_row,
                                    T0_row,
                                    t_grid,
                                    xgrid,
                                    K,
                                    g = function(z) 2 * atan(z) / pi,
                                    project = TRUE) {
  Phi <- create_cosine_basis(t_grid, K)
  Zg <- as.numeric(Phi %*% Xi_row)
  Ug <- g(Zg)
  L <- length(t_grid)
  G <- length(xgrid)
  out <- matrix(NA_real_, nrow = L, ncol = G)
  for (ell in seq_len(L)) {
    out[ell, ] <- scalar_mult_transport(Ug[ell], T0_row, xgrid, project = project)
  }
  out
}

imse_transports <- function(T_hat_list, T_true_list, t_grid, xgrid,
                            squared = FALSE,
                            return_pointwise = FALSE) {
  n <- length(T_hat_list)
  if (length(T_true_list) != n) stop("Length mismatch in imse_transports.")
  t_grid <- as.numeric(t_grid)
  xgrid  <- as.numeric(xgrid)
  L <- length(t_grid)
  if (L < 1) stop("imse_transports: t_grid must have positive length.")

  # Pointwise W1 error (or W1^2 if squared=TRUE): n x L
  err_mat <- matrix(NA_real_, nrow = n, ncol = L)

  for (i in seq_len(n)) {
    Th <- as.matrix(T_hat_list[[i]])
    Tt <- as.matrix(T_true_list[[i]])
    if (!all(dim(Th) == dim(Tt))) stop(sprintf("Subject %d: dim mismatch in T_hat vs T_true.", i))
    if (nrow(Th) != L) stop(sprintf("Subject %d: nrow(T_hat) must equal length(t_grid).", i))

    w1_t <- trapz_rows(xgrid, abs(Th - Tt))
    if (isTRUE(squared)) w1_t <- w1_t^2
    err_mat[i, ] <- w1_t
  }

  if (isTRUE(return_pointwise)) return(err_mat)

  # Integrated error over t for each subject
  trapz_rows(t_grid, err_mat)
}

# =========================
# End of file
# =========================
