#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

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

parse_args <- function(args) {
  out <- list()
  for (a in args) {
    if (!grepl("=", a, fixed = TRUE)) next
    key <- sub("=.*$", "", a)
    val <- sub("^[^=]*=", "", a)
    out[[key]] <- val
  }
  out
}

`%||%` <- function(x, y) if (!is.null(x)) x else y

monotone_quantile <- function(q) {
  q <- as.numeric(q)
  q <- cummax(q)
  q[1] <- 0
  q[length(q)] <- 1
  q
}

approx_inv_cdf <- function(qbar, u_grid, xgrid) {
  qbar <- monotone_quantile(qbar)
  # Remove duplicates in x to keep approx stable.
  keep <- !duplicated(qbar)
  x <- qbar[keep]
  y <- u_grid[keep]
  if (length(x) < 2) return(rep(0.5, length(xgrid)))
  ux <- approx(x = x, y = y, xout = xgrid, rule = 2, ties = "ordered")$y
  pmin(pmax(ux, 0), 1)
}

project_transport <- function(Tx) {
  Tx <- pmin(pmax(as.numeric(Tx), 0), 1)
  Tx <- cummax(Tx)
  Tx[1] <- 0
  Tx[length(Tx)] <- 1
  Tx
}

make_subject_transport_list <- function(df_roi_long, xgrid) {
  subjects <- sort(unique(df_roi_long$subject_id))
  T_list <- vector("list", length(subjects))
  t_list <- vector("list", length(subjects))
  obs_info <- vector("list", length(subjects))

  for (i in seq_along(subjects)) {
    sid <- subjects[i]
    dsub <- df_roi_long %>%
      filter(subject_id == sid) %>%
      select(subject_id, condition, run, t, x, T) %>%
      distinct() %>%
      arrange(t, x)

    obs <- dsub %>%
      select(condition, run, t) %>%
      distinct() %>%
      arrange(t)

    if (nrow(obs) == 0) {
      T_list[[i]] <- matrix(NA_real_, nrow = 0, ncol = length(xgrid))
      t_list[[i]] <- numeric(0)
      obs_info[[i]] <- obs
      next
    }

    T_mat <- matrix(NA_real_, nrow = nrow(obs), ncol = length(xgrid))
    for (j in seq_len(nrow(obs))) {
      tj <- obs$t[j]
      rowj <- dsub %>% filter(abs(t - tj) < 1e-12) %>% arrange(x)
      if (nrow(rowj) != length(xgrid)) {
        stop(sprintf("Subject %s time %.3f has %d x-points, expected %d.",
                     sid, tj, nrow(rowj), length(xgrid)))
      }
      T_mat[j, ] <- rowj$T
    }

    T_list[[i]] <- T_mat
    t_list[[i]] <- obs$t
    obs_info[[i]] <- obs
  }

  list(
    subjects = subjects,
    T_list = T_list,
    t_list = t_list,
    obs_info = obs_info
  )
}

main <- function() {
  kv <- parse_args(commandArgs(trailingOnly = TRUE))

  input_quant <- hcp_path(kv$input_quant %||% "data_processed/hcp_roi_quantiles_raw.csv")
  output_quant_norm <- hcp_path(kv$output_quant_norm %||% "data_processed/hcp_distribution_quantiles_normalized.csv")
  output_transports <- hcp_path(kv$output_transports %||% "data_processed/hcp_transports_long.csv")
  output_mean_quant <- hcp_path(kv$output_mean_quant %||% "data_processed/hcp_mean_quantiles.csv")
  output_objects <- hcp_path(kv$output_objects %||% "data_processed/hcp_transport_objects.rds")
  output_norm_stats <- hcp_path(kv$output_norm_stats %||% "data_processed/hcp_roi_norm_stats.csv")
  log_file <- hcp_path(kv$log_file %||% "logs/03_build_transports.log")
  xgrid_len <- as.integer(kv$xgrid_len %||% "201")

  if (!file.exists(input_quant)) {
    stop(sprintf("Input quantile file not found: %s", input_quant))
  }
  if (!is.finite(xgrid_len) || xgrid_len < 11) {
    stop("xgrid_len must be >= 11.")
  }

  for (path in c(log_file, output_quant_norm, output_transports, output_mean_quant, output_objects, output_norm_stats)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  }
  log_con <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)

  log_msg <- function(...) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
    writeLines(msg, con = log_con)
    cat(msg, "\n")
  }

  log_msg("Loading raw ROI quantiles from ", input_quant)
  qraw <- read_csv(input_quant, show_col_types = FALSE)
  required_cols <- c("subject_id", "condition", "run", "t", "roi", "p", "q_raw")
  if (!all(required_cols %in% names(qraw))) {
    stop("Input quantile file missing required columns.")
  }

  qraw <- qraw %>%
    mutate(
      subject_id = as.character(subject_id),
      condition = as.character(condition),
      run = as.character(run),
      roi = as.character(roi),
      t = as.numeric(t),
      p = as.numeric(p),
      q_raw = as.numeric(q_raw)
    ) %>%
    filter(is.finite(p), is.finite(q_raw), is.finite(t))

  norm_stats <- qraw %>%
    group_by(roi) %>%
    summarise(
      q_min = min(q_raw, na.rm = TRUE),
      q_max = max(q_raw, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(scale = ifelse(q_max > q_min, q_max - q_min, 1))

  qnorm <- qraw %>%
    left_join(norm_stats, by = "roi") %>%
    mutate(
      q_norm = (q_raw - q_min) / scale,
      q_norm = pmin(pmax(q_norm, 0), 1)
    ) %>%
    select(subject_id, condition, run, t, roi, p, q_raw, q_norm)

  write_csv(norm_stats, output_norm_stats)
  write_csv(qnorm, output_quant_norm)
  log_msg("Wrote normalization statistics: ", output_norm_stats)
  log_msg("Wrote normalized quantiles: ", output_quant_norm)

  u_grid <- sort(unique(qnorm$p))
  xgrid <- seq(0, 1, length.out = xgrid_len)
  rois <- sort(unique(qnorm$roi))
  log_msg("ROIs: ", paste(rois, collapse = ", "))
  log_msg("Number of quantile probabilities: ", length(u_grid), "; transport x-grid size: ", length(xgrid))

  mean_quant_rows <- list()
  transport_rows <- list()
  row_idx_mq <- 1L
  row_idx_tr <- 1L

  for (roi in rois) {
    d_roi <- qnorm %>% filter(roi == !!roi)
    keys <- d_roi %>%
      select(condition, run, t) %>%
      distinct() %>%
      arrange(t)

    for (kk in seq_len(nrow(keys))) {
      cond_k <- keys$condition[kk]
      run_k <- keys$run[kk]
      t_k <- keys$t[kk]

      d_key <- d_roi %>%
        filter(condition == cond_k, run == run_k, abs(t - t_k) < 1e-12) %>%
        select(subject_id, p, q_norm)

      wide <- d_key %>%
        pivot_wider(names_from = p, values_from = q_norm)
      p_cols <- setdiff(names(wide), "subject_id")
      p_num <- as.numeric(p_cols)
      ord <- order(p_num)
      p_cols <- p_cols[ord]
      p_num <- p_num[ord]

      q_mat <- as.matrix(wide[, p_cols, drop = FALSE])
      q_bar <- colMeans(q_mat, na.rm = TRUE)
      q_bar <- monotone_quantile(q_bar)
      u_of_x <- approx_inv_cdf(qbar = q_bar, u_grid = p_num, xgrid = xgrid)

      mean_quant_rows[[row_idx_mq]] <- data.frame(
        roi = roi,
        condition = cond_k,
        run = run_k,
        t = t_k,
        p = p_num,
        q_mean = q_bar
      )
      row_idx_mq <- row_idx_mq + 1L

      for (ss in seq_len(nrow(wide))) {
        sid <- wide$subject_id[ss]
        q_sub <- as.numeric(wide[ss, p_cols, drop = TRUE])
        q_sub <- monotone_quantile(q_sub)
        Tx <- approx(x = p_num, y = q_sub, xout = u_of_x, rule = 2, ties = "ordered")$y
        Tx <- project_transport(Tx)

        transport_rows[[row_idx_tr]] <- data.frame(
          subject_id = sid,
          roi = roi,
          condition = cond_k,
          run = run_k,
          t = t_k,
          x = xgrid,
          T = Tx
        )
        row_idx_tr <- row_idx_tr + 1L
      }
    }
  }

  mean_quant_df <- bind_rows(mean_quant_rows)
  transport_df <- bind_rows(transport_rows)

  write_csv(mean_quant_df, output_mean_quant)
  write_csv(transport_df, output_transports)
  log_msg("Wrote mean quantiles: ", output_mean_quant)
  log_msg("Wrote transport maps (long format): ", output_transports)

  transport_objects <- list()
  for (roi in rois) {
    d_roi <- transport_df %>% filter(roi == !!roi)
    built <- make_subject_transport_list(d_roi, xgrid = xgrid)
    transport_objects[[roi]] <- list(
      roi = roi,
      xgrid = xgrid,
      subjects = built$subjects,
      T_list = built$T_list,
      t_list = built$t_list,
      obs_info = built$obs_info
    )
  }

  saveRDS(
    list(
      transport_objects = transport_objects,
      mean_quantiles = mean_quant_df,
      norm_stats = norm_stats,
      u_grid = u_grid,
      xgrid = xgrid
    ),
    file = output_objects
  )
  log_msg("Saved transport objects RDS: ", output_objects)
  log_msg("Build transports completed.")
}

main()
