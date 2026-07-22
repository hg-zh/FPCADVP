# baseline_fpca_methods.R
# ------------------------------------------------------------
# Comparable baseline methods for distribution-/transport-valued
# processes (for simulation + real data):
#
# (1) LQD baseline (Petersen & Muller 2016b idea):
#     - For each observed distribution, compute log quantile density
#       y(u) = log Q'(u).
#     - Expand y(u) in a basis over u.
#     - For each basis coefficient, run longitudinal FPCA (PACE) over time.
#     - Reconstruct y_hat(t,u), then invert to quantile Q_hat(t,u).
#
# (2) Tangent/log-map Wasserstein FPCA baseline:
#     - Use log map at the reference (barycenter) distribution:
#       v(u) = T(u) - u  (since log_mu(nu) = T - id in 1D).
#     - Expand v(u) in a basis over u.
#     - Run longitudinal FPCA on each coefficient over time.
#     - Reconstruct v_hat(t,u); map back via exp: T_hat = id + v_hat,
#       and project to the transport space (monotone + bounds).
#
# Notes:
# - These are "blocked" main functions designed to be re-used in real data.
# - They accept list inputs: curves (quantiles or transports) and times.
# - Requires: fdapace, splines (recommended base package)
# ------------------------------------------------------------

options(stringsAsFactors = FALSE)

# =========================
# 0) Utilities (define only if not already loaded)
# =========================

if (!exists("clamp")) {
  clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)
}

if (!exists("project_transport")) {
  project_transport <- function(T, lower = 0, upper = 1) {
    T <- as.numeric(T)
    if (length(T) < 2) return(T)
    T <- clamp(T, lower, upper)
    T[1] <- lower
    T[length(T)] <- upper
    T <- cummax(T)
    T[1] <- lower
    T[length(T)] <- upper
    T
  }
}

if (!exists("trapz_vec")) {
  trapz_vec <- function(x, y) {
    x <- as.numeric(x); y <- as.numeric(y)
    if (length(x) != length(y)) stop("trapz_vec: length mismatch")
    if (length(x) < 2) return(0)
    dx <- diff(x)
    sum((y[-length(y)] + y[-1]) * dx / 2)
  }
}

# Row-wise trapezoid integration for Y (nrow x length(x))
if (!exists("trapz_rows")) {
  trapz_rows <- function(x, Y) {
    x <- as.numeric(x)
    Y <- as.matrix(Y)
    G <- length(x)
    if (ncol(Y) != G) stop("trapz_rows: grid mismatch")
    if (G < 2) return(rep(0, nrow(Y)))
    dx <- diff(x)
    mids <- (Y[, -G, drop = FALSE] + Y[, -1, drop = FALSE]) / 2
    mids <- sweep(mids, 2, dx, "*")
    rowSums(mids)
  }
}

# Cumulative trapezoid integration for each row of Y.
# Returns matrix nrow(Y) x length(x) with cumint at each x.
if (!exists("cumtrapz_rows")) {
  cumtrapz_rows <- function(x, Y) {
    x <- as.numeric(x)
    Y <- as.matrix(Y)
    G <- length(x)
    if (ncol(Y) != G) stop("cumtrapz_rows: grid mismatch")
    if (G < 2) return(matrix(0, nrow = nrow(Y), ncol = G))
    dx <- diff(x)  # length G-1
    incr <- (Y[, -G, drop = FALSE] + Y[, -1, drop = FALSE]) / 2
    incr <- sweep(incr, 2, dx, "*")
    out <- matrix(0, nrow = nrow(Y), ncol = G)
    for (k in 2:G) out[, k] <- out[, k - 1] + incr[, k - 1]
    out
  }
}

# =========================
# 1) Basis over u-grid (for representing functions in u)
# =========================

create_cosine_basis_01 <- function(grid, K) {
  grid <- as.numeric(grid)
  if (K < 1) stop("K must be >=1")
  Phi <- matrix(0, nrow = length(grid), ncol = K)
  Phi[, 1] <- 1
  if (K >= 2) {
    Phi[, 2:K] <- cos(2 * pi * outer(grid, 1:(K - 1)))
  }
  Phi
}

make_u_basis <- function(grid,
                         K,
                         type = c("bspline", "cosine"),
                         bs_degree = 3) {
  type <- match.arg(type)
  grid <- as.numeric(grid)

  if (type == "bspline") {
    if (!requireNamespace("splines", quietly = TRUE)) {
      stop("Package 'splines' is required for bspline basis (it is part of base R 'recommended').")
    }
    B <- splines::bs(grid, df = K, degree = bs_degree, intercept = TRUE)
    B <- as.matrix(B)
  } else {
    B <- create_cosine_basis_01(grid, K)
  }

  list(
    B = B,
    qrB = qr(B),
    type = type,
    K = K,
    grid = grid
  )
}

basis_coef_matrix <- function(curves_mat, qrB) {
  curves_mat <- as.matrix(curves_mat)
  Y <- t(curves_mat)  # G x Ni
  coefs <- qr.coef(qrB, Y)
  t(as.matrix(coefs))
}

# =========================
# 2) Scalar FPCA helper (PACE) + interpolation to a chosen t_grid
# =========================

fpca_predict_scalar <- function(Ly,
                                Lt,
                                t_grid = NULL,
                                methodXi = "CE",
                                fpca_optns = NULL,
                                save_fit = FALSE) {
  if (!requireNamespace("fdapace", quietly = TRUE)) {
    stop("Package 'fdapace' is required. Install it via install.packages('fdapace').")
  }

  if (is.null(fpca_optns)) fpca_optns <- list()
  fpca_optns$methodXi <- methodXi
  if (is.null(fpca_optns$methodMuCovEst)) fpca_optns$methodMuCovEst <- "smooth"

  fit <- fdapace::FPCA(Ly, Lt, optns = fpca_optns)

  workGrid <- as.numeric(fit$workGrid)
  mu <- as.numeric(fit$mu)
  traj <- fit$xiEst %*% t(fit$phi)
  traj <- sweep(traj, 2, mu, "+")

  # Interpolate to user grid if provided
  if (!is.null(t_grid)) {
    t_grid <- as.numeric(t_grid)
    n <- nrow(traj)
    traj2 <- matrix(NA_real_, nrow = n, ncol = length(t_grid))
    for (i in seq_len(n)) {
      traj2[i, ] <- approx(x = workGrid, y = traj[i, ], xout = t_grid, rule = 2)$y
    }
    traj <- traj2
    workGrid <- t_grid
  }

  out <- list(pred = traj, t_grid = workGrid)
  if (save_fit) out$fpca_fit <- fit
  out
}

# =========================
# 3) LQD transform + inverse
# =========================

# Central/forward/backward differences for derivative on [0,1]
quantile_derivative_fd <- function(Qmat, ugrid) {
  Qmat <- as.matrix(Qmat)
  ugrid <- as.numeric(ugrid)
  G <- length(ugrid)
  if (ncol(Qmat) != G) stop("quantile_derivative_fd: grid mismatch")
  if (G < 3) stop("Need at least 3 grid points for stable finite-difference derivative")

  dQ <- matrix(NA_real_, nrow = nrow(Qmat), ncol = G)

  # boundaries: forward/backward
  dQ[, 1] <- (Qmat[, 2] - Qmat[, 1]) / (ugrid[2] - ugrid[1])
  dQ[, G] <- (Qmat[, G] - Qmat[, G - 1]) / (ugrid[G] - ugrid[G - 1])

  # central differences
  denom <- (ugrid[3:G] - ugrid[1:(G - 2)])
  dQ[, 2:(G - 1)] <- (Qmat[, 3:G] - Qmat[, 1:(G - 2)]) / rep(denom, each = nrow(Qmat))

  dQ
}

lqd_transform_quantiles <- function(Qmat,
                                    ugrid,
                                    eps_deriv = 1e-8,
                                    smooth_deriv = FALSE,
                                    spar = NULL) {
  Qmat <- as.matrix(Qmat)
  ugrid <- as.numeric(ugrid)
  G <- length(ugrid)
  if (ncol(Qmat) != G) stop("lqd_transform_quantiles: grid mismatch")

  if (!smooth_deriv) {
    dQ <- quantile_derivative_fd(Qmat, ugrid)
  } else {
    # smooth.spline per curve, then derivative
    dQ <- matrix(NA_real_, nrow = nrow(Qmat), ncol = G)
    for (j in seq_len(nrow(Qmat))) {
      sp <- if (is.null(spar)) {
        stats::smooth.spline(ugrid, Qmat[j, ])
      } else {
        stats::smooth.spline(ugrid, Qmat[j, ], spar = spar)
      }
      dQ[j, ] <- stats::predict(sp, ugrid, deriv = 1)$y
    }
  }

  dQ <- pmax(dQ, eps_deriv)
  log(dQ)
}

lqd_inverse_logqd <- function(logqd_mat,
                              ugrid,
                              lower = 0,
                              upper = 1,
                              project = TRUE) {
  logqd_mat <- as.matrix(logqd_mat)
  ugrid <- as.numeric(ugrid)
  G <- length(ugrid)
  if (ncol(logqd_mat) != G) stop("lqd_inverse_logqd: grid mismatch")

  qd <- exp(logqd_mat)
  cum <- cumtrapz_rows(ugrid, qd)  # each row integrates from 0
  denom <- cum[, G]
  denom[!is.finite(denom) | denom <= 0] <- 1
  Q <- cum / denom

  # map to [lower, upper]
  Q <- lower + (upper - lower) * Q

  # enforce exact endpoints
  Q[, 1] <- lower
  Q[, G] <- upper

  if (project) {
    for (j in seq_len(nrow(Q))) {
      Q[j, ] <- project_transport(Q[j, ], lower = lower, upper = upper)
    }
  }

  Q
}

# =========================
# 4) Main baseline: LQD + FPCA over time
# =========================

# Q_list: list of Ni x G matrices of quantile functions Q_ij(u) evaluated on ugrid.
# In your simulation, you can pass T_hat because T_hat is the empirical quantile.
fit_LQD_fpca <- function(Q_list,
                         t_list,
                         ugrid,
                         Ku = 10,
                         basis_type = c("bspline", "cosine"),
                         bs_degree = 3,
                         t_grid = seq(0, 1, length.out = 51),
                         methodXi = "CE",
                         fpca_optns = NULL,
                         smooth_deriv = FALSE,
                         spar = NULL,
                         eps_deriv = 1e-8,
                         lower = NULL,
                         upper = NULL,
                         project = TRUE,
                         save_fpca_fits = FALSE,
                         verbose = FALSE) {

  if (!is.list(Q_list) || !is.list(t_list)) stop("Q_list and t_list must be lists.")
  n <- length(Q_list)
  if (length(t_list) != n) stop("Length mismatch: Q_list and t_list must match.")

  ugrid <- as.numeric(ugrid)
  G <- length(ugrid)

  # Reorder each subject by time
  for (i in seq_len(n)) {
    ti <- as.numeric(t_list[[i]])
    Qi <- as.matrix(Q_list[[i]])
    if (ncol(Qi) != G) stop(sprintf("Subject %d: ncol(Q_i) != length(ugrid)", i))
    if (length(ti) != nrow(Qi)) stop(sprintf("Subject %d: length(t_i) != nrow(Q_i)", i))
    o <- order(ti)
    t_list[[i]] <- ti[o]
    Q_list[[i]] <- Qi[o, , drop = FALSE]
  }

  # Default support endpoints (if not provided): use global min/max of observed endpoints.
  if (is.null(lower)) {
    lower <- min(unlist(lapply(Q_list, function(Qi) Qi[, 1])), na.rm = TRUE)
  }
  if (is.null(upper)) {
    upper <- max(unlist(lapply(Q_list, function(Qi) Qi[, ncol(Qi)])), na.rm = TRUE)
  }

  # Basis in u
  basis_type <- match.arg(basis_type)
  ub <- make_u_basis(ugrid, Ku, type = basis_type, bs_degree = bs_degree)
  B <- ub$B
  qrB <- ub$qrB

  # Transform each observed quantile curve to log quantile density
  coef_list <- vector("list", n)
  for (i in seq_len(n)) {
    if (verbose && (i %% 50 == 0)) cat(sprintf("  [LQD] transform subject %d/%d\n", i, n))
    logqd_i <- lqd_transform_quantiles(
      Qmat = Q_list[[i]],
      ugrid = ugrid,
      eps_deriv = eps_deriv,
      smooth_deriv = smooth_deriv,
      spar = spar
    )
    coef_list[[i]] <- basis_coef_matrix(logqd_i, qrB)  # Ni x Ku
  }

  # FPCA per coefficient
  L <- length(t_grid)
  coef_hat <- array(NA_real_, dim = c(n, L, Ku))
  fpca_fits <- if (save_fpca_fits) vector("list", Ku) else NULL

  for (k in seq_len(Ku)) {
    if (verbose) cat(sprintf("  [LQD] FPCA coef %d/%d\n", k, Ku))
    Ly_k <- lapply(seq_len(n), function(i) as.numeric(coef_list[[i]][, k]))
    pred_k <- fpca_predict_scalar(
      Ly = Ly_k,
      Lt = t_list,
      t_grid = t_grid,
      methodXi = methodXi,
      fpca_optns = fpca_optns,
      save_fit = save_fpca_fits
    )
    coef_hat[, , k] <- pred_k$pred
    if (save_fpca_fits) fpca_fits[[k]] <- pred_k$fpca_fit
  }

  # Reconstruct logqd and invert to quantiles
  Q_hat_list <- vector("list", n)
  for (i in seq_len(n)) {
    coef_mat <- matrix(coef_hat[i, , ], nrow = L, ncol = Ku)
    logqd_hat <- coef_mat %*% t(B)  # L x G
    Q_hat <- lqd_inverse_logqd(logqd_hat, ugrid, lower = lower, upper = upper, project = project)
    Q_hat_list[[i]] <- Q_hat
  }

  structure(list(
    method = "LQD_FPCA",
    ugrid = ugrid,
    t_grid = as.numeric(t_grid),
    Ku = Ku,
    basis_type = basis_type,
    basis = B,
    coef_hat = coef_hat,
    Q_hat = Q_hat_list,
    fpca_fits = fpca_fits,
    call = match.call()
  ), class = "lqd_fpca_fit")
}

# =========================
# 5) Main baseline: Tangent/log-map Wasserstein FPCA
# =========================

# T_list: list of Ni x G matrices of transport maps T_ij(x) evaluated on xgrid.
# This assumes transports are from the reference (barycenter) to each subject,
# so log-map is v = T - id.
fit_tangent_wass_fpca <- function(T_list,
                                  t_list,
                                  xgrid,
                                  Ku = 10,
                                  basis_type = c("bspline", "cosine"),
                                  bs_degree = 3,
                                  t_grid = seq(0, 1, length.out = 51),
                                  methodXi = "CE",
                                  fpca_optns = NULL,
                                  lower = NULL,
                                  upper = NULL,
                                  project = TRUE,
                                  save_fpca_fits = FALSE,
                                  verbose = FALSE) {

  if (!is.list(T_list) || !is.list(t_list)) stop("T_list and t_list must be lists.")
  n <- length(T_list)
  if (length(t_list) != n) stop("Length mismatch: T_list and t_list must match.")

  xgrid <- as.numeric(xgrid)
  G <- length(xgrid)

  # Reorder each subject by time
  for (i in seq_len(n)) {
    ti <- as.numeric(t_list[[i]])
    Ti <- as.matrix(T_list[[i]])
    if (ncol(Ti) != G) stop(sprintf("Subject %d: ncol(T_i) != length(xgrid)", i))
    if (length(ti) != nrow(Ti)) stop(sprintf("Subject %d: length(t_i) != nrow(T_i)", i))
    o <- order(ti)
    t_list[[i]] <- ti[o]
    T_list[[i]] <- Ti[o, , drop = FALSE]
  }

  # Default transport bounds from xgrid (if not provided)
  if (is.null(lower)) lower <- min(xgrid)
  if (is.null(upper)) upper <- max(xgrid)

  # Basis in x
  basis_type <- match.arg(basis_type)
  ub <- make_u_basis(xgrid, Ku, type = basis_type, bs_degree = bs_degree)
  B <- ub$B
  qrB <- ub$qrB

  # Project observed transports to tangent vectors v = T - id
  coef_list <- vector("list", n)
  for (i in seq_len(n)) {
    if (verbose && (i %% 50 == 0)) cat(sprintf("  [Tangent] prep subject %d/%d\n", i, n))
    Vi <- sweep(as.matrix(T_list[[i]]), 2, xgrid, "-")  # Ni x G
    coef_list[[i]] <- basis_coef_matrix(Vi, qrB)        # Ni x Ku
  }

  # FPCA per coefficient
  L <- length(t_grid)
  coef_hat <- array(NA_real_, dim = c(n, L, Ku))
  fpca_fits <- if (save_fpca_fits) vector("list", Ku) else NULL

  for (k in seq_len(Ku)) {
    if (verbose) cat(sprintf("  [Tangent] FPCA coef %d/%d\n", k, Ku))
    Ly_k <- lapply(seq_len(n), function(i) as.numeric(coef_list[[i]][, k]))
    pred_k <- fpca_predict_scalar(
      Ly = Ly_k,
      Lt = t_list,
      t_grid = t_grid,
      methodXi = methodXi,
      fpca_optns = fpca_optns,
      save_fit = save_fpca_fits
    )
    coef_hat[, , k] <- pred_k$pred
    if (save_fpca_fits) fpca_fits[[k]] <- pred_k$fpca_fit
  }

  # Reconstruct v_hat and map back to transports T_hat = id + v_hat
  T_hat_list <- vector("list", n)
  for (i in seq_len(n)) {
    coef_mat <- matrix(coef_hat[i, , ], nrow = L, ncol = Ku)
    v_hat <- coef_mat %*% t(B)  # L x G
    T_hat <- sweep(v_hat, 2, xgrid, "+")
    if (project) {
      for (ell in seq_len(nrow(T_hat))) {
        T_hat[ell, ] <- project_transport(T_hat[ell, ], lower = lower, upper = upper)
      }
    }
    T_hat_list[[i]] <- T_hat
  }

  structure(list(
    method = "Tangent_Wasserstein_FPCA",
    xgrid = xgrid,
    t_grid = as.numeric(t_grid),
    Ku = Ku,
    basis_type = basis_type,
    basis = B,
    coef_hat = coef_hat,
    T_hat = T_hat_list,
    fpca_fits = fpca_fits,
    call = match.call()
  ), class = "tangent_fpca_fit")
}

# =========================
# End of file
# =========================
