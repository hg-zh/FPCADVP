# transport_fpca_modular.R
# ------------------------------------------------------------
# Modular implementation for "transport-based FPCA" for
# distribution-valued (transport-valued) processes.
#
# Designed so the SAME main function can be used for
# - simulations (where T_ij are estimated from samples), and
# - real data (where T_ij are estimated upstream).
#
# The structure follows the estimation steps in Sections 3.1–3.3
# of your submission.  (Dense vs sparse is controlled by `latent`.)
# ------------------------------------------------------------

# =========================
# 0) Small helpers
# =========================

clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

# Trapezoidal integral for a vector y evaluated on grid x
trapz_vec <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) != length(y)) stop("trapz_vec: length mismatch")
  if (length(x) < 2) return(0)
  dx <- diff(x)
  sum((y[-length(y)] + y[-1]) * dx / 2)
}

# Row-wise trapezoidal integral for a matrix Y: n x length(x)
trapz_rows <- function(x, Y) {
  x <- as.numeric(x)
  Y <- as.matrix(Y)
  G <- length(x)
  if (ncol(Y) != G) stop("trapz_rows: ncol(Y) must equal length(x)")
  if (G < 2) return(rep(0, nrow(Y)))
  dx <- diff(x)  # length G-1
  mids <- (Y[, -G, drop = FALSE] + Y[, -1, drop = FALSE]) / 2
  mids <- sweep(mids, 2, dx, "*")
  rowSums(mids)
}

# Enforce [0,1] bounds and monotonicity (nondecreasing).
# (This matches your Get_order + boundary clipping, but in vectorized form.)
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

# Fast inverse of a monotone transport map T evaluated on xgrid.
# Returns T^{-1} evaluated on the SAME xgrid.
# (Much faster than pointwise uniroot.)
inverse_transport <- function(T, xgrid) {
  T <- as.numeric(T)
  xgrid <- as.numeric(xgrid)
  if (length(T) != length(xgrid)) stop("inverse_transport: length mismatch")
  Tm <- project_transport(T, lower = min(xgrid), upper = max(xgrid))

  # Handle ties in x=T by collapsing to unique points
  idx <- !duplicated(Tm)
  x_u <- Tm[idx]
  y_u <- xgrid[idx]
  if (length(x_u) < 2) return(xgrid)  # degenerate fallback
  approx(x = x_u, y = y_u, xout = xgrid, rule = 2, ties = "ordered")$y
}

# sign(T) = sign( ∫ (T(u) - u) du ) on the grid
transport_sign_rows <- function(Tmat, xgrid, tol = 1e-12) {
  xgrid <- as.numeric(xgrid)
  Tmat <- as.matrix(Tmat)
  if (ncol(Tmat) != length(xgrid)) stop("transport_sign_rows: grid mismatch")
  disp <- sweep(Tmat, 2, xgrid, "-")
  ints <- trapz_rows(xgrid, disp)
  s <- sign(ints)
  # To avoid empty sets from numerical near-zeros, treat |int|<tol as +1 by default
  s[abs(ints) < tol] <- 1
  s
}

# ||T||_1 = ∫ |T(u) - u| du on the grid
transport_norm1_rows <- function(Tmat, xgrid) {
  xgrid <- as.numeric(xgrid)
  Tmat <- as.matrix(Tmat)
  disp <- abs(sweep(Tmat, 2, xgrid, "-"))
  trapz_rows(xgrid, disp)
}

# W1 distance between two transports on the grid
w1_distance_vec <- function(T1, T2, xgrid) {
  xgrid <- as.numeric(xgrid)
  trapz_vec(xgrid, abs(as.numeric(T1) - as.numeric(T2)))
}

# Rescale a transport so that ∫|T(u)-u|du equals `target_norm`
# (Your equation (8)-style rescaling; we project back to [0,1] afterwards.)
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

# α ⊙ T0 for α in [-1,1], using your scalar multiplication definition.
# If you want to avoid negative α, you can pass T0_minus (inverse-class baseline)
# and use alpha=abs(U) (see main function below).
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
  # T_i: Ni x G matrix of observed/estimated transports for subject i
  T_i <- as.matrix(T_i)
  signs <- transport_sign_rows(T_i, xgrid, tol = tol_sign)
  I_pos <- which(signs > 0)
  I_neg <- which(signs < 0)

  # Raw Fréchet integrals (just means in L2) in each sign group
  T_plus_raw  <- if (length(I_pos) > 0) colMeans(T_i[I_pos, , drop = FALSE]) else rep(NA_real_, length(xgrid))
  T_minus_raw <- if (length(I_neg) > 0) colMeans(T_i[I_neg, , drop = FALSE]) else rep(NA_real_, length(xgrid))

  # Rescale to have norm kappa
  T_plus  <- if (!anyNA(T_plus_raw))  rescale_transport_to_norm(T_plus_raw,  xgrid, kappa, project = project) else rep(NA_real_, length(xgrid))
  T_minus <- if (!anyNA(T_minus_raw)) rescale_transport_to_norm(T_minus_raw, xgrid, kappa, project = project) else rep(NA_real_, length(xgrid))

  # If one side is missing, approximate it by inverting the other side
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

estimate_U_hat <- function(T_i, signs, xgrid, kappa, eps_u = 1e-6) {
  norms <- transport_norm1_rows(T_i, xgrid)
  u <- as.numeric(signs) * norms / kappa
  clamp(u, -1 + eps_u, 1 - eps_u)
}

# =========================
# 2) Main method: transport FPCA
# =========================
# Inputs:
# - T_list: list length n; element i is Ni x G matrix (rows correspond to times t_list[[i]])
# - t_list: list length n; element i is length-Ni vector of times in [0,1]
# - xgrid : transport grid (u-grid), e.g. seq(0,1,length=128)
#
# Key knobs:
# - latent=TRUE  -> sparse/longitudinal method: run FPCA on Z = g^{-1}(U) and then map back via g
# - latent=FALSE -> dense-style method: run FPCA directly on U (may need clamping)
#
fit_transport_fpca <- function(T_list,
                               t_list,
                               xgrid,
                               kappa = NULL,
                               latent = TRUE,
                               g = function(z) 2 * atan(z) / pi,
                               g_inv = function(u) tan(pi * u / 2),
                               fpca_optns = NULL,
                               methodXi = NULL,
                               eps_u = 1e-6,
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

  # Reorder each subject's data so times are increasing (safer for FPCA)
  for (i in seq_len(n)) {
    ti <- as.numeric(t_list[[i]])
    Ti <- as.matrix(T_list[[i]])
    if (ncol(Ti) != G) stop(sprintf("Subject %d: ncol(T_i) must equal length(xgrid).", i))
    if (length(ti) != nrow(Ti)) stop(sprintf("Subject %d: length(t_i) must equal nrow(T_i).", i))
    o <- order(ti)
    t_list[[i]] <- ti[o]
    T_list[[i]] <- Ti[o, , drop = FALSE]
  }

  # Choose kappa (fixed baseline norm) if not supplied.
  # A safe default is the maximum observed transport norm so that |U_hat| <= 1.
  if (is.null(kappa)) {
    all_norms <- unlist(lapply(T_list, function(Ti) transport_norm1_rows(Ti, xgrid)))
    kappa <- max(all_norms, na.rm = TRUE)
    if (!is.finite(kappa) || kappa <= 0) stop("Could not infer a valid kappa from the data.")
    if (verbose) message(sprintf("kappa not supplied; using max observed norm: kappa = %.4f", kappa))
  }

  # 1) Estimate subject-specific baselines T_{i0} (positive and negative classes)
  T0_plus  <- matrix(NA_real_, nrow = n, ncol = G)
  T0_minus <- matrix(NA_real_, nrow = n, ncol = G)
  signs_list <- vector("list", n)

  for (i in seq_len(n)) {
    base_i <- estimate_subject_baselines(T_list[[i]], xgrid, kappa, project = project, tol_sign = tol_sign)
    T0_plus[i, ]  <- base_i$T_plus
    T0_minus[i, ] <- base_i$T_minus
    signs_list[[i]] <- base_i$signs
  }

  # 2) Estimate U_hat(t_ij) from norms/signs
  U_hat_list <- vector("list", n)
  for (i in seq_len(n)) {
    U_hat_list[[i]] <- estimate_U_hat(T_list[[i]], signs_list[[i]], xgrid, kappa, eps_u = eps_u)
  }

  # 3) Construct scalar observations for FPCA:
  #    - latent=TRUE : Z_hat = g^{-1}(U_hat), then FPCA on Z_hat
  #    - latent=FALSE: FPCA on U_hat directly
  if (latent) {
    Ly <- lapply(U_hat_list, g_inv)
  } else {
    Ly <- U_hat_list
  }
  Lt <- t_list

  # 4) Fit FPCA (PACE)
  if (is.null(methodXi)) methodXi <- if (latent) "CE" else "IN"
  if (is.null(fpca_optns)) fpca_optns <- list()
  fpca_optns$methodXi <- methodXi
  if (is.null(fpca_optns$methodMuCovEst)) fpca_optns$methodMuCovEst <- "smooth"

  fpca_fit <- fdapace::FPCA(Ly, Lt, optns = fpca_optns)

  # 5) Reconstruct scalar trajectories on workGrid
  workGrid <- as.numeric(fpca_fit$workGrid)
  L <- length(workGrid)
  mu <- as.numeric(fpca_fit$mu)
  if (length(mu) != L) {
    # Some fdapace versions return mu as matrix; try to coerce
    mu <- as.numeric(mu)
    if (length(mu) != L) stop("fdapace::FPCA returned mu with unexpected length.")
  }

  # xiEst: n x K, phi: L x K  => xiEst %*% t(phi) is n x L
  traj_hat <- fpca_fit$xiEst %*% t(fpca_fit$phi)
  traj_hat <- sweep(traj_hat, 2, mu, "+")

  if (latent) {
    Z_hat_grid <- traj_hat
    U_hat_grid <- g(Z_hat_grid)  # guaranteed in (-1,1)
  } else {
    Z_hat_grid <- NULL
    U_hat_grid <- clamp(traj_hat, -1 + eps_u, 1 - eps_u)
  }

  # 6) Reconstruct transport trajectories on workGrid
  T_hat_list <- vector("list", n)
  for (i in seq_len(n)) {
    Ti_hat <- matrix(NA_real_, nrow = L, ncol = G)
    for (ell in seq_len(L)) {
      u_val <- U_hat_grid[i, ell]
      if (u_val >= 0) {
        Ti_hat[ell, ] <- xgrid + u_val * (T0_plus[i, ] - xgrid)
      } else {
        # Use baseline in the negative equivalence class and the positive magnitude |u_val|
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
    eps_u = eps_u,
    t_list = t_list,
    T_obs = T_list,
    signs_obs = signs_list,
    T0_plus = T0_plus,
    T0_minus = T0_minus,
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
# 3) Simulation utilities
# =========================

# Cosine basis used in your simulation section:
# phi_1(t)=1, phi_k(t)=cos(2*(k-1)*pi*t) for k>=2.
create_cosine_basis <- function(t, K) {
  t <- as.numeric(t)
  if (K < 1) stop("K must be >= 1")

  # For exact reproducibility with the original simulation code (sim_sparse_illu.R),
  # we use fdapace::CreateBasis(type="cos") when it is available.
  if (requireNamespace("fdapace", quietly = TRUE)) {
    # CreateBasis is exported by fdapace (used directly in the legacy script)
    return(fdapace::CreateBasis(K, t, type = "cos"))
  }

  # Fallback: plain cosine basis (matches the manuscript definition)
  Phi <- matrix(0, nrow = length(t), ncol = K)
  Phi[, 1] <- 1
  if (K >= 2) {
    Phi[, 2:K] <- cos(2 * pi * outer(t, 1:(K - 1)))
  }
  Phi
}


# Sample m draws from the distribution with quantile (transport) function T(u) on grid xgrid.
sample_from_transport <- function(T, xgrid, m) {
  u <- runif(m, 0, 1)
  approx(x = xgrid, y = as.numeric(T), xout = u, rule = 2)$y
}

# Estimate a transport map (quantile function) from a numeric vector of samples.
estimate_transport_from_samples <- function(samples, probs, enforce_bounds = TRUE) {
  q <- as.numeric(stats::quantile(samples, probs = probs, type = 8, names = FALSE, na.rm = TRUE))
  if (enforce_bounds) {
    q[1] <- 0
    q[length(q)] <- 1
    q <- project_transport(q)
  }
  q
}

# Estimate all T_ij for one subject from an Ni x m matrix of samples
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

# Generate one synthetic dataset (simulation replication)
simulate_transport_process <- function(n,
                                       N,
                                       m,
                                       xgrid = seq(0, 1, length.out = 128),
                                       K = 50,
                                       t_design = c("random", "fixed"),
                                       g = function(z) 2 * atan(z) / pi,
                                       beta_a = c(3, 4),
                                       beta_b = c(1, 2),
                                       seed = NULL,
                                       project = TRUE) {

  t_design <- match.arg(t_design)
  if (!is.null(seed)) set.seed(seed)

  xgrid <- as.numeric(xgrid)
  G <- length(xgrid)

  # Observation times
  t_list <- vector("list", n)
  if (t_design == "random") {
    for (i in seq_len(n)) t_list[[i]] <- sort(runif(N, 0, 1))
  } else {
    # equispaced on (0,1)
    tg <- seq(0, 1, length.out = N + 2)[2:(N + 1)]
    for (i in seq_len(n)) t_list[[i]] <- tg
  }

  # Latent Gaussian coefficients.
# IMPORTANT: To reproduce the original Figure 2 simulation (sim_sparse_illu.R),
# we must (i) generate coefficients row-wise (i then k) to match the RNG stream,
# and (ii) use sd = 1/k^1.5 (as in the original code).
#
# This is equivalent to the nested loop:
#   for (i in 1:n) for (k in 1:K) Xi[i,k] <- (1/k^1.5) * rnorm(1)
Xi_raw <- matrix(rnorm(n * K, mean = 0, sd = 1), nrow = n, ncol = K, byrow = TRUE)
Xi <- sweep(Xi_raw, 2, (1:K)^(1.5), "/")

  # Baseline transports T_i0 as Beta quantiles.
# IMPORTANT: The original sim_sparse_illu.R draws (a_i, b_i) *within* the subject loop
# (a_1, b_1, a_2, b_2, ...). We replicate this RNG order for exact figure reproduction.
a <- numeric(n)
b <- numeric(n)
T0_raw <- matrix(NA_real_, nrow = n, ncol = G)
for (i in seq_len(n)) {
  a[i] <- runif(1, beta_a[1], beta_a[2])
  b[i] <- runif(1, beta_b[1], beta_b[2])
  T0_raw[i, ] <- stats::qbeta(xgrid, a[i], b[i])
}

  norms0 <- transport_norm1_rows(T0_raw, xgrid)
  kappa_true <- mean(norms0)
  T0 <- matrix(NA_real_, nrow = n, ncol = G)
  for (i in seq_len(n)) T0[i, ] <- rescale_transport_to_norm(T0_raw[i, ], xgrid, kappa_true, project = project)

  # Generate Z_i(t_ij), U_i(t_ij)=g(Z_i(t_ij)), true transports, samples, and estimated transports
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
      # True transport at observed time
      Ti_true[j, ] <- scalar_mult_transport(Ui[j], T0[i, ], xgrid, project = project)
      # Sample from the induced distribution
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

# Compute true transports on an arbitrary time grid for the simulated model
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

# =========================
# 4) Evaluation helper (IMSE)
# =========================

# IMSE_i = ∫ W1(T_hat_i(t), T_true_i(t)) dt on a common grid t_grid,
# approximated by trapezoidal integration.
imse_transports <- function(T_hat_list, T_true_list, t_grid, xgrid) {
  n <- length(T_hat_list)
  if (length(T_true_list) != n) stop("Length mismatch in imse_transports.")
  out <- numeric(n)
  for (i in seq_len(n)) {
    Th <- as.matrix(T_hat_list[[i]])
    Tt <- as.matrix(T_true_list[[i]])
    if (!all(dim(Th) == dim(Tt))) stop(sprintf("Subject %d: dim mismatch in T_hat vs T_true.", i))
    # W1 at each time point
    errs <- apply(abs(Th - Tt), 1, function(row_abs) trapz_vec(xgrid, row_abs))
    out[i] <- trapz_vec(t_grid, errs)
  }
  out
}

# =========================
# End of file
# =========================
