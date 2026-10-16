#!/usr/bin/env Rscript

# RStudio-friendly generator for outputs/figures/hcp_eigentransport_modes.pdf.
# Run run_hcp_eigentransport_fit.R first if the fit object does not exist.

suppressPackageStartupMessages({
  library(ggplot2)
})

# ==============================
# User-editable configuration
# ==============================

t_rest <- 0.40
t_motor <- 0.925

fit_rds <- NULL
out_pdf <- NULL



# Paths default to this application directory, regardless of the working directory.
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file)) {
  script_file <- sub("^--file=", "", script_file[1])
} else {
  script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(script_file) && requireNamespace("rstudioapi", quietly = TRUE)) {
    script_file <- rstudioapi::getSourceEditorContext()$path
  }
}
if (is.null(script_file) || !nzchar(script_file)) stop("Run this file with Rscript or Source.")
script_dir <- dirname(normalizePath(gsub("~+~", " ", script_file, fixed = TRUE), mustWork = TRUE))
source(file.path(script_dir, "hcp_paths.R"))
kv <- hcp_parse_args(commandArgs(trailingOnly = TRUE))
project_root <- kv$root %||% script_dir
hcp_path <- function(path) hcp_resolve_path(path, project_root)


if (is.null(fit_rds)) {
  fit_rds <- hcp_path(kv$fit_rds %||% "outputs/hcp_transport_fpca_fits.rds")
}
if (is.null(out_pdf)) {
  out_pdf <- hcp_path(kv$out_pdf %||% "outputs/figures/hcp_eigentransport_modes.pdf")
}

log_file <- file.path(project_root, "logs", "make_hcp_eigentransport_modes.log")
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

log_con <- file(log_file, open = "wt")
on.exit(close(log_con), add = TRUE)

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  writeLines(msg, con = log_con)
  cat(msg, "\n")
}

project_transport <- function(Tx) {
  Tx <- pmin(pmax(as.numeric(Tx), 0), 1)
  Tx <- cummax(Tx)
  Tx[1] <- 0
  Tx[length(Tx)] <- 1
  Tx
}

if (!file.exists(fit_rds)) {
  stop("Fit object not found at ", fit_rds, ". Run run_hcp_eigentransport_fit.R first.")
}

obj <- readRDS(fit_rds)
fit_by_roi <- obj$fit_by_roi
if (is.null(fit_by_roi) || length(fit_by_roi) == 0) {
  stop("No fitted ROI objects found in ", fit_rds)
}

log_msg("Building eigentransport mode figure from ", fit_rds)

mode_rows <- list()
idx <- 1L

for (roi in names(fit_by_roi)) {
  fit <- fit_by_roi[[roi]]$fit
  xgrid <- fit$xgrid
  work_grid <- fit$workGrid
  mu <- as.numeric(fit$fpca_fit$mu)
  phi <- as.matrix(fit$fpca_fit$phi)

  if (ncol(phi) < 1) {
    log_msg("Skipping ROI ", roi, " because no eigenfunctions were fitted.")
    next
  }

  lambda1 <- as.numeric(fit$fpca_fit$lambda[1])
  if (!is.finite(lambda1) || lambda1 < 0) {
    log_msg("Skipping ROI ", roi, " because lambda1 is invalid.")
    next
  }

  sd1 <- sqrt(lambda1)
  gfun <- fit$g
  Tplus_mean <- colMeans(fit$T0_plus, na.rm = TRUE)
  Tminus_mean <- colMeans(fit$T0_minus, na.rm = TRUE)

  time_lookup <- c(Rest = t_rest, Motor = t_motor)
  for (time_label in names(time_lookup)) {
    tt <- time_lookup[[time_label]]
    mu_t <- approx(x = work_grid, y = mu, xout = tt, rule = 2)$y
    phi_t <- approx(x = work_grid, y = phi[, 1], xout = tt, rule = 2)$y

    z_vals <- c(mu_t - 2 * sd1 * phi_t, mu_t, mu_t + 2 * sd1 * phi_t)
    z_labels <- c("PC1 - 2SD", "Mean", "PC1 + 2SD")

    for (ii in seq_along(z_vals)) {
      u <- as.numeric(gfun(z_vals[ii]))
      if (u >= 0) {
        Tx <- xgrid + u * (Tplus_mean - xgrid)
      } else {
        Tx <- xgrid + (-u) * (Tminus_mean - xgrid)
      }
      Tx <- project_transport(Tx)

      mode_rows[[idx]] <- data.frame(
        roi = roi,
        time_label = time_label,
        mode = z_labels[ii],
        x = xgrid,
        transport_minus_id = Tx - xgrid
      )
      idx <- idx + 1L
    }
  }
}

mode_df <- do.call(rbind, mode_rows)
if (is.null(mode_df) || nrow(mode_df) == 0) {
  stop("No eigentransport curves were generated.")
}

p_mode <- ggplot(mode_df, aes(x = x, y = transport_minus_id, color = mode, linetype = mode)) +
  geom_line(linewidth = 0.8) +
  facet_grid(time_label ~ roi) +
  labs(
    x = "x",
    y = "T(x) - x"
  ) +
  scale_x_continuous(
    labels = function(x) ifelse(x %in% c(0, 1), as.character(x), sprintf("%.2f", x))
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  guides(
    color = guide_legend(nrow = 1),
    linetype = guide_legend(nrow = 1)
  )

p_mode

ggsave(out_pdf, p_mode, width = 800/72, height = 450/72)
log_msg("Saved figure to ", out_pdf)
