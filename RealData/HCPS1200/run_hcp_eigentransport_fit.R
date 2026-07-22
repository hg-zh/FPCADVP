#!/usr/bin/env Rscript

# Minimal RStudio-friendly fit script for the HCP eigentransport analysis.
# Dependencies:
#   1. data_processed/hcp_transport_objects.rds
#   2. proposed_method/transport_fpca_modular.R
# Output:
#   outputs/hcp_transport_fpca_fits.rds

# ==============================
# User-editable configuration
# ==============================

input_transport_rds <- NULL
method_script <- NULL
output_fit_rds <- NULL


project_root <- "/Users/hangzhou/Documents/research/Hans/random object/JRSSB_RR1/Final_Code/RealData/HCPS1200"
setwd(project_root)

if (is.null(input_transport_rds)) {
  input_transport_rds <- file.path(project_root, "data_processed", "hcp_transport_objects.rds")
}
if (is.null(method_script)) {
  method_script <- file.path(project_root, "transport_fpca_modular.R")
}
if (is.null(output_fit_rds)) {
  output_fit_rds <- file.path(project_root, "outputs", "hcp_transport_fpca_fits.rds")
}

log_file <- file.path(project_root, "logs", "run_hcp_eigentransport_fit.log")
dir.create(dirname(output_fit_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_con <- file(log_file, open = "wt")
on.exit(close(log_con), add = TRUE)

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  writeLines(msg, con = log_con)
  cat(msg, "\n")
}

if (!file.exists(input_transport_rds)) {
  stop("Missing transport object file: ", input_transport_rds)
}
if (!file.exists(method_script)) {
  stop("Missing proposed method script: ", method_script)
}

log_msg("Loading transport objects from ", input_transport_rds)
transport_obj <- readRDS(input_transport_rds)
transport_objects <- transport_obj$transport_objects

if (is.null(transport_objects) || length(transport_objects) == 0) {
  stop("No transport objects found in ", input_transport_rds)
}

source(method_script)

fit_one_roi <- function(obj, roi_name) {
  keep <- which(vapply(obj$T_list, nrow, integer(1)) > 0)
  if (length(keep) == 0) {
    log_msg("Skipping ROI ", roi_name, " because it has no valid subject trajectories.")
    return(NULL)
  }

  log_msg("Fitting transport-FPCA for ROI: ", roi_name)
  fit <- fit_transport_fpca(
    T_list = obj$T_list[keep],
    t_list = obj$t_list[keep],
    xgrid = obj$xgrid,
    latent = TRUE,
    methodXi = "CE",
    u_cap = 0.95,
    Z_cap = 8,
    fpca_optns = list(maxK = 3, FVEthreshold = 0.999),
    verbose = FALSE
  )

  list(
    fit = fit,
    subjects = obj$subjects[keep],
    obs_info = obj$obs_info[keep]
  )
}

fit_by_roi <- list()
for (roi in names(transport_objects)) {
  fit_res <- fit_one_roi(transport_objects[[roi]], roi)
  if (!is.null(fit_res)) {
    fit_by_roi[[roi]] <- fit_res
  }
}

if (length(fit_by_roi) == 0) {
  stop("No ROI fits were produced.")
}

saveRDS(
  list(
    fit_by_roi = fit_by_roi,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    input_transport_rds = normalizePath(input_transport_rds, winslash = "/", mustWork = TRUE),
    method_script = normalizePath(method_script, winslash = "/", mustWork = TRUE)
  ),
  file = output_fit_rds
)

log_msg("Saved fitted object to ", output_fit_rds)
log_msg("Next step: source('scripts/make_hcp_eigentransport_modes.R') to rebuild the figure.")
