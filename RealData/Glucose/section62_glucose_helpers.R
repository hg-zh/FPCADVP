`%||%` <- function(x, y) if (is.null(x)) y else x

glucose_clamp <- function(x, lo, hi) {
  pmin(pmax(x, lo), hi)
}

glucose_enforce_transport <- function(x, lo, hi) {
  x <- as.numeric(x)
  x <- glucose_clamp(x, lo, hi)
  x[1] <- lo
  x[length(x)] <- hi
  x <- cummax(x)
  x[1] <- lo
  x[length(x)] <- hi
  x
}

glucose_quantile_to_cdf <- function(prob_grid, quantile_values, support_grid) {
  quantile_values <- as.numeric(quantile_values)
  idx <- !duplicated(quantile_values)
  stats::approx(
    x = quantile_values[idx],
    y = prob_grid[idx],
    xout = support_grid,
    rule = 2,
    ties = "ordered"
  )$y
}

glucose_transport_from_mean_cdf <- function(mean_cdf, prob_grid, subject_quantile) {
  stats::approx(
    x = prob_grid,
    y = subject_quantile,
    xout = mean_cdf,
    rule = 2,
    ties = "ordered"
  )$y
}

glucose_load_data <- function(data_path, day_window = 10) {
  data_raw <- read.csv(data_path, stringsAsFactors = FALSE)
  # Match the legacy Section 6.2 preprocessing, which parsed timestamps in the
  # local R session timezone before converting to study-day offsets.
  data_raw$time <- as.POSIXct(data_raw$time)
  data_raw$gl <- suppressWarnings(as.numeric(data_raw$gl))
  data_raw <- stats::na.omit(data_raw)

  subject_ids <- unique(data_raw$id)
  data_raw$time_res <- 0

  for (i in seq_along(subject_ids)) {
    idx <- data_raw$id == subject_ids[i]
    dates <- as.Date(data_raw$time[idx])
    data_raw$time_res[idx] <- as.numeric(dates - dates[1])
  }

  keep_days <- vector("list", length(subject_ids))
  names(keep_days) <- subject_ids
  for (i in seq_along(subject_ids)) {
    candidate_days <- sort(unique(data_raw$time_res[data_raw$id == subject_ids[i]]))
    keep_days[[i]] <- candidate_days[candidate_days - candidate_days[1] <= day_window]
  }

  list(
    data = data_raw,
    subject_ids = subject_ids,
    keep_days = keep_days
  )
}

glucose_make_split <- function(keep_days, seed, holdout_days = 2, min_days_for_holdout = 6) {
  set.seed(seed)

  n_subjects <- length(keep_days)
  train_days <- vector("list", n_subjects)
  test_days <- vector("list", n_subjects)

  for (i in seq_along(keep_days)) {
    days_i <- keep_days[[i]]
    if (length(days_i) >= min_days_for_holdout) {
      held_out <- sort(sample(days_i, holdout_days))
      train_days[[i]] <- setdiff(days_i, held_out)
      test_days[[i]] <- held_out
    } else {
      train_days[[i]] <- days_i
      test_days[[i]] <- numeric(0)
    }
  }

  list(train_days = train_days, test_days = test_days)
}

glucose_make_day_samples <- function(data_raw, subject_ids, keep_days) {
  out <- vector("list", length(subject_ids))
  names(out) <- subject_ids

  for (i in seq_along(subject_ids)) {
    entries <- vector("list", length(keep_days[[i]]))
    names(entries) <- as.character(keep_days[[i]])
    for (j in seq_along(keep_days[[i]])) {
      day_j <- keep_days[[i]][j]
      entries[[j]] <- data_raw$gl[data_raw$id == subject_ids[i] & data_raw$time_res == day_j]
    }
    out[[i]] <- entries
  }

  out
}

glucose_subject_quantiles <- function(day_samples,
                                      support_grid,
                                      prob_grid,
                                      spline_spar = 0.5) {
  out <- vector("list", length(day_samples))
  names(out) <- names(day_samples)

  lo <- min(support_grid)
  hi <- max(support_grid)

  for (i in seq_along(day_samples)) {
    entries <- day_samples[[i]]
    q_mat <- matrix(
      NA_real_,
      nrow = length(entries),
      ncol = length(prob_grid),
      dimnames = list(names(entries), NULL)
    )

    for (j in seq_along(entries)) {
      q0 <- as.numeric(stats::quantile(entries[[j]], probs = prob_grid, names = FALSE, na.rm = TRUE))
      q0 <- glucose_clamp(q0, lo, hi)
      q0[1] <- lo
      q0[length(q0)] <- hi

      q_smooth <- stats::smooth.spline(prob_grid, q0, spar = spline_spar)$y
      q_smooth <- glucose_enforce_transport(q_smooth, lo, hi)
      q_mat[j, ] <- q_smooth
    }

    out[[i]] <- q_mat
  }

  out
}

glucose_fit_mean <- function(day_samples,
                             subject_ids,
                             train_days,
                             support_grid,
                             time_grid,
                             bw_reg = 4,
                             kernel_reg = "epan",
                             den_lower_threshold = 0.01,
                             min_points_per_day = 30) {
  xin <- unlist(train_days)
  while (length(unique(xin)) != length(xin)) {
    xin <- xin + stats::runif(length(xin), 0, 1e-6)
  }

  yin <- list()
  for (i in seq_along(subject_ids)) {
    for (day_i in train_days[[i]]) {
      yin[[length(yin) + 1]] <- day_samples[[i]][[as.character(day_i)]]
    }
  }

  keep <- which(vapply(yin, length, integer(1)) > min_points_per_day)
  xin <- xin[keep]
  yin <- yin[keep]

  mean_fit <- frechet::LocDenReg(
    xin = xin,
    yin = yin,
    xout = time_grid,
    optns = list(
      bwReg = bw_reg,
      kernelReg = kernel_reg,
      denLowerThreshold = den_lower_threshold,
      dSup = support_grid
    )
  )

  mean_quantiles <- mean_fit$qout
  mean_quantiles <- pmin(pmax(mean_quantiles, min(support_grid)), max(support_grid))
  mean_quantiles[, 1] <- min(support_grid)
  mean_quantiles[, ncol(mean_quantiles)] <- max(support_grid)

  mean_cdf <- t(apply(mean_quantiles, 1, function(q_row) {
    cdf_row <- glucose_quantile_to_cdf(mean_fit$qSup, q_row, support_grid)
    cdf_row[1] <- 0
    cdf_row[length(cdf_row)] <- 1
    cummax(cdf_row)
  }))

  list(
    mean_fit = mean_fit,
    mean_quantiles = mean_quantiles,
    mean_cdf = mean_cdf
  )
}

glucose_make_transports <- function(subject_quantiles,
                                    keep_days,
                                    mean_cdf,
                                    time_grid,
                                    support_grid) {
  support_lo <- min(support_grid)
  support_hi <- max(support_grid)
  support_width <- support_hi - support_lo
  prob_grid <- seq(0, 1, length.out = ncol(subject_quantiles[[1]]))

  orig <- vector("list", length(subject_quantiles))
  scaled <- vector("list", length(subject_quantiles))
  names(orig) <- names(subject_quantiles)
  names(scaled) <- names(subject_quantiles)

  for (i in seq_along(subject_quantiles)) {
    days_i <- keep_days[[i]]
    q_mat <- subject_quantiles[[i]]
    T_orig <- matrix(
      NA_real_,
      nrow = length(days_i),
      ncol = length(support_grid),
      dimnames = list(as.character(days_i), NULL)
    )

    for (j in seq_along(days_i)) {
      time_index <- which(abs(time_grid - days_i[j]) < 1e-10)
      transport_j <- glucose_transport_from_mean_cdf(mean_cdf[time_index, ], prob_grid, q_mat[as.character(days_i[j]), ])
      T_orig[j, ] <- glucose_enforce_transport(transport_j, support_lo, support_hi)
    }

    orig[[i]] <- T_orig
    scaled[[i]] <- (T_orig - support_lo) / support_width
  }

  list(orig = orig, scaled = scaled)
}

glucose_pick_rows <- function(T_list, days_list) {
  out <- vector("list", length(T_list))
  names(out) <- names(T_list)

  for (i in seq_along(T_list)) {
    if (length(days_list[[i]]) == 0) {
      out[[i]] <- matrix(numeric(0), nrow = 0, ncol = ncol(T_list[[i]]))
    } else {
      out[[i]] <- T_list[[i]][as.character(days_list[[i]]), , drop = FALSE]
    }
  }

  out
}

glucose_predict_at_times <- function(fit, target_times, support_grid) {
  support_lo <- min(support_grid)
  support_hi <- max(support_grid)
  support_width <- support_hi - support_lo

  out <- vector("list", length(fit$T_hat))

  for (i in seq_along(fit$T_hat)) {
    times_i <- target_times[[i]]
    if (length(times_i) == 0) {
      out[[i]] <- matrix(numeric(0), nrow = 0, ncol = length(fit$xgrid))
      next
    }

    pred_scaled <- matrix(
      NA_real_,
      nrow = length(times_i),
      ncol = length(fit$xgrid),
      dimnames = list(as.character(times_i), NULL)
    )

    for (g in seq_along(fit$xgrid)) {
      pred_scaled[, g] <- stats::approx(
        x = fit$workGrid,
        y = fit$T_hat[[i]][, g],
        xout = times_i,
        rule = 2,
        ties = "ordered"
      )$y
    }

    for (r in seq_len(nrow(pred_scaled))) {
      pred_scaled[r, ] <- project_transport(pred_scaled[r, ])
    }

    out[[i]] <- pred_scaled * support_width + support_lo
    rownames(out[[i]]) <- as.character(times_i)
  }

  out
}

glucose_quantile_predictions_to_transports <- function(Q_hat_list,
                                                       mean_cdf,
                                                       time_grid,
                                                       support_grid,
                                                       prob_grid = NULL) {
  if (is.null(prob_grid)) {
    prob_grid <- seq(0, 1, length.out = ncol(Q_hat_list[[1]]))
  }

  support_lo <- min(support_grid)
  support_hi <- max(support_grid)
  out <- vector("list", length(Q_hat_list))

  for (i in seq_along(Q_hat_list)) {
    Q_hat_i <- as.matrix(Q_hat_list[[i]])
    T_hat_i <- matrix(
      NA_real_,
      nrow = nrow(Q_hat_i),
      ncol = length(support_grid),
      dimnames = list(as.character(time_grid), NULL)
    )

    for (ell in seq_len(nrow(Q_hat_i))) {
      T_hat_i[ell, ] <- glucose_enforce_transport(
        glucose_transport_from_mean_cdf(
          mean_cdf = mean_cdf[ell, ],
          prob_grid = prob_grid,
          subject_quantile = Q_hat_i[ell, ]
        ),
        support_lo,
        support_hi
      )
    }

    out[[i]] <- T_hat_i
  }

  out
}

glucose_select_prediction_times <- function(curve_list, time_grid, target_times) {
  out <- vector("list", length(curve_list))

  for (i in seq_along(curve_list)) {
    times_i <- target_times[[i]]
    curve_i <- as.matrix(curve_list[[i]])

    if (length(times_i) == 0) {
      out[[i]] <- matrix(numeric(0), nrow = 0, ncol = ncol(curve_i))
      next
    }

    idx <- match(times_i, time_grid)
    if (all(!is.na(idx))) {
      out[[i]] <- curve_i[idx, , drop = FALSE]
      rownames(out[[i]]) <- as.character(times_i)
      next
    }

    pred_i <- matrix(
      NA_real_,
      nrow = length(times_i),
      ncol = ncol(curve_i),
      dimnames = list(as.character(times_i), NULL)
    )
    for (g in seq_len(ncol(curve_i))) {
      pred_i[, g] <- stats::approx(
        x = time_grid,
        y = curve_i[, g],
        xout = times_i,
        rule = 2,
        ties = "ordered"
      )$y
    }
    out[[i]] <- pred_i
  }

  out
}

glucose_holdout_error_table <- function(pred_orig, obs_orig, subject_ids, test_days, support_grid) {
  rows <- list()

  for (i in seq_along(subject_ids)) {
    if (length(test_days[[i]]) == 0) {
      next
    }

    pointwise_errors <- numeric(length(test_days[[i]]))
    for (j in seq_along(test_days[[i]])) {
      day_j <- test_days[[i]][j]
      pointwise_errors[j] <- w1_distance_vec(
        pred_orig[[i]][as.character(day_j), ],
        obs_orig[[i]][as.character(day_j), ],
        support_grid
      )
    }

    obs_signs <- transport_sign_rows(obs_orig[[i]][as.character(test_days[[i]]), , drop = FALSE], support_grid)
    pred_signs <- transport_sign_rows(pred_orig[[i]][as.character(test_days[[i]]), , drop = FALSE], support_grid)

    rows[[length(rows) + 1]] <- data.frame(
      subject_index = i,
      subject_id = subject_ids[i],
      n_holdout_days = length(test_days[[i]]),
      holdout_days = paste(test_days[[i]], collapse = ","),
      holdout_errors = paste(sprintf("%.6f", pointwise_errors), collapse = ","),
      mean_holdout_w1 = mean(pointwise_errors),
      observed_signs = paste(obs_signs, collapse = ","),
      predicted_signs = paste(pred_signs, collapse = ","),
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  out[order(out$mean_holdout_w1), ]
}

glucose_default_panel_titles <- function(selected_summary, multiline = FALSE) {
  panel_titles <- character(0)

  for (k in seq_len(nrow(selected_summary))) {
    days_k <- as.numeric(strsplit(selected_summary$holdout_days[k], ",", fixed = TRUE)[[1]])
    subject_index <- selected_summary$subject_index[k]
    subject_id <- selected_summary$subject_id[k]

    for (day_j in days_k) {
      if (isTRUE(multiline)) {
        panel_titles <- c(
          panel_titles,
          paste0("Subject ", subject_index, "\n(", subject_id, "), t = ", day_j)
        )
      } else {
        panel_titles <- c(
          panel_titles,
          paste0("Subject ", subject_index, " (", subject_id, "), t = ", day_j)
        )
      }
    }
  }

  panel_titles
}

glucose_make_plot_df <- function(selected_summary,
                                 pred_orig,
                                 obs_orig,
                                 support_grid,
                                 panel_titles = NULL) {
  rows <- list()
  default_panel_titles <- glucose_default_panel_titles(selected_summary, multiline = FALSE)
  panel_counter <- 0

  for (k in seq_len(nrow(selected_summary))) {
    days_k <- as.numeric(strsplit(selected_summary$holdout_days[k], ",", fixed = TRUE)[[1]])
    subject_index <- selected_summary$subject_index[k]

    for (day_j in days_k) {
      panel_counter <- panel_counter + 1
      panel_label <- default_panel_titles[panel_counter]

      rows[[length(rows) + 1]] <- data.frame(
        subject_rank = k,
        subject_index = subject_index,
        subject_id = selected_summary$subject_id[k],
        holdout_day = day_j,
        panel = panel_label,
        glucose = support_grid,
        displacement = obs_orig[[subject_index]][as.character(day_j), ] - support_grid,
        curve_type = "Observed",
        stringsAsFactors = FALSE
      )

      rows[[length(rows) + 1]] <- data.frame(
        subject_rank = k,
        subject_index = subject_index,
        subject_id = selected_summary$subject_id[k],
        holdout_day = day_j,
        panel = panel_label,
        glucose = support_grid,
        displacement = pred_orig[[subject_index]][as.character(day_j), ] - support_grid,
        curve_type = "Predicted",
        stringsAsFactors = FALSE
      )
    }
  }

  plot_df <- do.call(rbind, rows)

  final_panel_titles <- default_panel_titles
  if (!is.null(panel_titles)) {
    if (length(panel_titles) != length(default_panel_titles)) {
      stop(sprintf(
        "panel_titles must have length %s, but got %s.",
        length(default_panel_titles),
        length(panel_titles)
      ))
    }
    final_panel_titles <- as.character(panel_titles)
  }

  panel_map <- stats::setNames(final_panel_titles, default_panel_titles)
  plot_df$panel <- unname(panel_map[plot_df$panel])
  plot_df$panel <- factor(plot_df$panel, levels = final_panel_titles)
  plot_df$curve_type <- factor(plot_df$curve_type, levels = c("Observed", "Predicted"))
  plot_df
}

glucose_make_method_comparison_df <- function(selected_summary,
                                              observed_holdouts,
                                              method_holdouts,
                                              support_grid,
                                              panel_titles = NULL) {
  if (is.null(names(method_holdouts)) || any(names(method_holdouts) == "")) {
    stop("method_holdouts must be a named list, e.g. list(Proposed = ..., LQD = ..., Tangent = ...).")
  }

  rows <- list()
  default_panel_titles <- glucose_default_panel_titles(selected_summary, multiline = FALSE)
  panel_counter <- 0

  for (k in seq_len(nrow(selected_summary))) {
    days_k <- as.numeric(strsplit(selected_summary$holdout_days[k], ",", fixed = TRUE)[[1]])
    subject_index <- selected_summary$subject_index[k]

    for (day_j in days_k) {
      panel_counter <- panel_counter + 1
      panel_label <- default_panel_titles[panel_counter]

      for (method_name in names(method_holdouts)) {
        rows[[length(rows) + 1]] <- data.frame(
          subject_rank = k,
          subject_index = subject_index,
          subject_id = selected_summary$subject_id[k],
          holdout_day = day_j,
          panel = panel_label,
          glucose = support_grid,
          displacement = method_holdouts[[method_name]][[subject_index]][as.character(day_j), ] - support_grid,
          curve_label = method_name,
          draw_order = 1L,
          stringsAsFactors = FALSE
        )
      }

      rows[[length(rows) + 1]] <- data.frame(
        subject_rank = k,
        subject_index = subject_index,
        subject_id = selected_summary$subject_id[k],
        holdout_day = day_j,
        panel = panel_label,
        glucose = support_grid,
        displacement = observed_holdouts[[subject_index]][as.character(day_j), ] - support_grid,
        curve_label = "Observed",
        draw_order = 2L,
        stringsAsFactors = FALSE
      )
    }
  }

  plot_df <- do.call(rbind, rows)

  final_panel_titles <- default_panel_titles
  if (!is.null(panel_titles)) {
    if (length(panel_titles) != length(default_panel_titles)) {
      stop(sprintf(
        "panel_titles must have length %s, but got %s.",
        length(default_panel_titles),
        length(panel_titles)
      ))
    }
    final_panel_titles <- as.character(panel_titles)
  }

  panel_map <- stats::setNames(final_panel_titles, default_panel_titles)
  plot_df$panel <- unname(panel_map[plot_df$panel])
  plot_df$panel <- factor(plot_df$panel, levels = final_panel_titles)
  plot_df$curve_label <- factor(
    plot_df$curve_label,
    levels = c("Observed", names(method_holdouts))
  )
  plot_df <- plot_df[order(plot_df$panel, plot_df$draw_order), ]
  rownames(plot_df) <- NULL
  plot_df
}

glucose_plot_selected_subjects <- function(selected_summary,
                                           pred_orig,
                                           obs_orig,
                                           support_grid,
                                           out_file,
                                           panel_titles = NULL,
                                           line_width = 0.5,
                                           zero_line_width = 0.3,
                                           axis_title_size = 10,
                                           axis_text_size = 9,
                                           strip_title_size = 10,
                                           legend_text_size = 9,
                                           legend_position = "bottom",
                                           legend_key_width = 1.2,
                                           legend_key_width_units = "cm",
                                           show_strip_text = TRUE,
                                           show_strip_background = TRUE,
                                           strip_background_fill = "white",
                                           strip_background_color = "black",
                                           strip_background_linewidth = 0.5,
                                           save_width = 13,
                                           save_height = 3.8,
                                           save_units = "in") {
  plot_df <- glucose_make_plot_df(
    selected_summary = selected_summary,
    pred_orig = pred_orig,
    obs_orig = obs_orig,
    support_grid = support_grid,
    panel_titles = panel_titles
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = glucose, y = displacement, linetype = curve_type)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = zero_line_width) +
    ggplot2::geom_line(color = "black", linewidth = line_width) +
    ggplot2::facet_wrap(~panel, nrow = 1) +
    ggplot2::scale_linetype_manual(values = c("Observed" = "solid", "Predicted" = "dashed")) +
    ggplot2::labs(
      x = "Glucose level (mg/dL)",
      y = "T(x)-x",
      linetype = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      legend.position = legend_position,
      axis.title = ggplot2::element_text(size = axis_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      strip.text = if (isTRUE(show_strip_text)) {
        ggplot2::element_text(size = strip_title_size)
      } else {
        ggplot2::element_blank()
      },
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.key.width = grid::unit(legend_key_width, legend_key_width_units),
      strip.background = if (isTRUE(show_strip_background)) {
        ggplot2::element_rect(
          fill = strip_background_fill,
          color = strip_background_color,
          linewidth = strip_background_linewidth
        )
      } else {
        ggplot2::element_blank()
      },
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(out_file, p, width = save_width, height = save_height, units = save_units)
  invisible(p)
}

glucose_plot_method_comparison <- function(selected_summary,
                                           observed_holdouts,
                                           method_holdouts,
                                           support_grid,
                                           out_file,
                                           panel_titles = NULL,
                                           method_colors = c(
                                             Observed = "black",
                                             Proposed = "#1f77b4",
                                             LQD = "#d95f02",
                                             Tangent = "#1b9e77"
                                           ),
                                           method_linetypes = c(
                                             Observed = "solid",
                                             Proposed = "dashed",
                                             LQD = "dotdash",
                                             Tangent = "longdash"
                                           ),
                                           line_width = 0.8,
                                           zero_line_width = 0.3,
                                           axis_title_size = 10,
                                           axis_text_size = 9,
                                           strip_title_size = 10,
                                           legend_text_size = 9,
                                           legend_position = "bottom",
                                           legend_key_width = 1.2,
                                           legend_key_width_units = "cm",
                                           show_strip_text = TRUE,
                                           show_strip_background = TRUE,
                                           strip_background_fill = "white",
                                           strip_background_color = "black",
                                           strip_background_linewidth = 0.5,
                                           save_width = 13,
                                           save_height = 3.8,
                                           save_units = "in") {
  plot_df <- glucose_make_method_comparison_df(
    selected_summary = selected_summary,
    observed_holdouts = observed_holdouts,
    method_holdouts = method_holdouts,
    support_grid = support_grid,
    panel_titles = panel_titles
  )

  legend_levels <- levels(plot_df$curve_label)
  if (!all(legend_levels %in% names(method_colors))) {
    stop("method_colors must be a named vector covering Observed and all methods.")
  }
  if (!all(legend_levels %in% names(method_linetypes))) {
    stop("method_linetypes must be a named vector covering Observed and all methods.")
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = glucose, y = displacement, color = curve_label, linetype = curve_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = zero_line_width) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::facet_wrap(~panel, nrow = 1) +
    ggplot2::scale_color_manual(values = method_colors[legend_levels]) +
    ggplot2::scale_linetype_manual(values = method_linetypes[legend_levels]) +
    ggplot2::labs(
      x = "Glucose level (mg/dL)",
      y = "Transport minus identity",
      color = NULL,
      linetype = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      legend.position = legend_position,
      axis.title = ggplot2::element_text(size = axis_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      strip.text = if (isTRUE(show_strip_text)) {
        ggplot2::element_text(size = strip_title_size)
      } else {
        ggplot2::element_blank()
      },
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.key.width = grid::unit(legend_key_width, legend_key_width_units),
      strip.background = if (isTRUE(show_strip_background)) {
        ggplot2::element_rect(
          fill = strip_background_fill,
          color = strip_background_color,
          linewidth = strip_background_linewidth
        )
      } else {
        ggplot2::element_blank()
      },
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(out_file, p, width = save_width, height = save_height, units = save_units)
  invisible(p)
}

glucose_write_report <- function(report_file,
                                 seed,
                                 selected_summary,
                                 kappa_original,
                                 kappa_scaled,
                                 fit_mode) {
  error_definition <- paste(
    "For subject i with held-out days t_1,...,t_m, the reconstruction error is",
    "mean_holdout_w1(i) = (1/m) * sum_l int_{50}^{160} |T_hat_i(t_l, x) - T_obs_i(t_l, x)| dx,",
    "computed numerically on the common 201-point support grid."
  )

  lines <- c(
    paste("Seed:", seed),
    paste("Fit mode:", fit_mode),
    paste("kappa_original:", sprintf("%.6f", kappa_original)),
    paste("kappa_scaled:", sprintf("%.10f", kappa_scaled)),
    paste("Error definition:", error_definition),
    ""
  )

  for (k in seq_len(nrow(selected_summary))) {
    lines <- c(
      lines,
      paste("Rank", k, ":"),
      paste("  subject_index =", selected_summary$subject_index[k]),
      paste("  subject_id =", selected_summary$subject_id[k]),
      paste("  holdout_days =", selected_summary$holdout_days[k]),
      paste("  holdout_errors =", selected_summary$holdout_errors[k]),
      paste("  mean_holdout_w1 =", sprintf("%.6f", selected_summary$mean_holdout_w1[k])),
      paste("  observed_signs =", selected_summary$observed_signs[k]),
      paste("  predicted_signs =", selected_summary$predicted_signs[k]),
      ""
    )
  }

  writeLines(lines, report_file)
}

glucose_run_section62 <- function(data_path,
                                  seed = 3,
                                  output_dir,
                                  day_window = 10,
                                  holdout_days = 2,
                                  min_days_for_holdout = 6,
                                  support_grid = seq(50, 160, length.out = 201),
                                  time_grid = seq(0, 10, by = 0.25),
                                  prob_grid = seq(0, 1, length.out = 201),
                                  bw_reg = 4,
                                  kernel_reg = "epan",
                                  den_lower_threshold = 0.01,
                                  min_points_per_day = 30,
                                  spline_spar = 0.5,
                                  kappa_original = 3200,
                                  latent = TRUE,
                                  u_cap = 0.95,
                                  Z_cap = 8,
                                  top_n = 2) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  support_lo <- min(support_grid)
  support_hi <- max(support_grid)
  support_width <- support_hi - support_lo
  xgrid_scaled <- (support_grid - support_lo) / support_width
  kappa_scaled <- kappa_original / (support_width^2)

  loaded <- glucose_load_data(data_path = data_path, day_window = day_window)
  split <- glucose_make_split(
    keep_days = loaded$keep_days,
    seed = seed,
    holdout_days = holdout_days,
    min_days_for_holdout = min_days_for_holdout
  )
  day_samples <- glucose_make_day_samples(loaded$data, loaded$subject_ids, loaded$keep_days)
  subject_quantiles <- glucose_subject_quantiles(day_samples, support_grid, prob_grid, spline_spar = spline_spar)
  mean_fit <- glucose_fit_mean(
    day_samples = day_samples,
    subject_ids = loaded$subject_ids,
    train_days = split$train_days,
    support_grid = support_grid,
    time_grid = time_grid,
    bw_reg = bw_reg,
    kernel_reg = kernel_reg,
    den_lower_threshold = den_lower_threshold,
    min_points_per_day = min_points_per_day
  )
  transports <- glucose_make_transports(
    subject_quantiles = subject_quantiles,
    keep_days = loaded$keep_days,
    mean_cdf = mean_fit$mean_cdf,
    time_grid = time_grid,
    support_grid = support_grid
  )

  fit <- fit_transport_fpca(
    T_list = glucose_pick_rows(transports$scaled, split$train_days),
    t_list = split$train_days,
    xgrid = xgrid_scaled,
    kappa = kappa_scaled,
    latent = latent,
    fpca_optns = list(
      nRegGrid = length(time_grid),
      outPercent = c(0, 1),
      error = FALSE
    ),
    methodXi = "CE",
    u_cap = u_cap,
    Z_cap = if (latent) Z_cap else NULL,
    project = TRUE
  )

  pred_orig <- glucose_predict_at_times(fit, split$test_days, support_grid)
  obs_orig <- glucose_pick_rows(transports$orig, split$test_days)
  error_table <- glucose_holdout_error_table(pred_orig, obs_orig, loaded$subject_ids, split$test_days, support_grid)
  selected_summary <- utils::head(error_table, top_n)
  selected_summary$rank <- seq_len(nrow(selected_summary))
  selected_summary <- selected_summary[, c(
    "rank",
    "subject_index",
    "subject_id",
    "n_holdout_days",
    "holdout_days",
    "holdout_errors",
    "mean_holdout_w1",
    "observed_signs",
    "predicted_signs"
  )]

  figure_file <- file.path(output_dir, sprintf("ind_gluco_seed_%s.pdf", seed))
  summary_file <- file.path(output_dir, sprintf("subject_errors_seed_%s.csv", seed))
  selected_file <- file.path(output_dir, sprintf("selected_subjects_seed_%s.csv", seed))
  report_file <- file.path(output_dir, sprintf("selection_report_seed_%s.txt", seed))
  results_file <- file.path(output_dir, sprintf("section62_results_seed_%s.rds", seed))

  glucose_plot_selected_subjects(selected_summary, pred_orig, obs_orig, support_grid, figure_file)
  utils::write.csv(error_table, summary_file, row.names = FALSE)
  utils::write.csv(selected_summary, selected_file, row.names = FALSE)
  glucose_write_report(
    report_file = report_file,
    seed = seed,
    selected_summary = selected_summary,
    kappa_original = kappa_original,
    kappa_scaled = kappa_scaled,
    fit_mode = if (latent) "latent_gaussian" else "direct_U"
  )

  saveRDS(
    list(
      seed = seed,
      config = list(
        day_window = day_window,
        holdout_days = holdout_days,
        min_days_for_holdout = min_days_for_holdout,
        bw_reg = bw_reg,
        kernel_reg = kernel_reg,
        den_lower_threshold = den_lower_threshold,
        min_points_per_day = min_points_per_day,
        spline_spar = spline_spar,
        kappa_original = kappa_original,
        kappa_scaled = kappa_scaled,
        latent = latent,
        u_cap = u_cap,
        Z_cap = if (latent) Z_cap else NA_real_
      ),
      subject_ids = loaded$subject_ids,
      keep_days = loaded$keep_days,
      split = split,
      mean_fit = mean_fit,
      fit = fit,
      error_table = error_table,
      selected_summary = selected_summary,
      observed_holdouts = obs_orig,
      predicted_holdouts = pred_orig,
      files = list(
        figure = figure_file,
        all_errors = summary_file,
        selected = selected_file,
        report = report_file
      )
    ),
    results_file
  )

  list(
    seed = seed,
    figure_file = figure_file,
    summary_file = summary_file,
    selected_file = selected_file,
    report_file = report_file,
    results_file = results_file,
    selected_summary = selected_summary,
    error_table = error_table,
    fit = fit
  )
}

glucose_classify_sign_string <- function(sign_string) {
  vals <- as.integer(strsplit(sign_string, ",", fixed = TRUE)[[1]])
  if (all(vals > 0)) return("positive")
  if (all(vals < 0)) return("negative")
  "mixed"
}

glucose_seed_illustration_table <- function(data_path,
                                            seeds,
                                            day_window = 10,
                                            holdout_days = 2,
                                            min_days_for_holdout = 6,
                                            support_grid = seq(50, 160, length.out = 201),
                                            time_grid = seq(0, 10, by = 0.25),
                                            prob_grid = seq(0, 1, length.out = 201),
                                            bw_reg = 4,
                                            kernel_reg = "epan",
                                            den_lower_threshold = 0.01,
                                            min_points_per_day = 30,
                                            spline_spar = 0.5,
                                            kappa_original = 3200,
                                            latent = TRUE,
                                            u_cap = 0.95,
                                            Z_cap = 8,
                                            top_n = 2,
                                            progress = FALSE) {
  loaded <- glucose_load_data(data_path = data_path, day_window = day_window)
  day_samples <- glucose_make_day_samples(loaded$data, loaded$subject_ids, loaded$keep_days)
  subject_quantiles <- glucose_subject_quantiles(
    day_samples = day_samples,
    support_grid = support_grid,
    prob_grid = prob_grid,
    spline_spar = spline_spar
  )

  support_width <- max(support_grid) - min(support_grid)
  xgrid_scaled <- (support_grid - min(support_grid)) / support_width
  kappa_scaled <- kappa_original / (support_width^2)

  rows <- list()

  for (seed in seeds) {
    if (isTRUE(progress)) {
      message(sprintf("Scoring seed %s", seed))
    }

    split <- glucose_make_split(
      keep_days = loaded$keep_days,
      seed = seed,
      holdout_days = holdout_days,
      min_days_for_holdout = min_days_for_holdout
    )

    mean_fit <- glucose_fit_mean(
      day_samples = day_samples,
      subject_ids = loaded$subject_ids,
      train_days = split$train_days,
      support_grid = support_grid,
      time_grid = time_grid,
      bw_reg = bw_reg,
      kernel_reg = kernel_reg,
      den_lower_threshold = den_lower_threshold,
      min_points_per_day = min_points_per_day
    )

    transports <- glucose_make_transports(
      subject_quantiles = subject_quantiles,
      keep_days = loaded$keep_days,
      mean_cdf = mean_fit$mean_cdf,
      time_grid = time_grid,
      support_grid = support_grid
    )

    fit <- fit_transport_fpca(
      T_list = glucose_pick_rows(transports$scaled, split$train_days),
      t_list = split$train_days,
      xgrid = xgrid_scaled,
      kappa = kappa_scaled,
      latent = latent,
      fpca_optns = list(
        nRegGrid = length(time_grid),
        outPercent = c(0, 1),
        error = FALSE
      ),
      methodXi = "CE",
      u_cap = u_cap,
      Z_cap = if (latent) Z_cap else NULL,
      project = TRUE
    )

    pred_orig <- glucose_predict_at_times(fit, split$test_days, support_grid)
    obs_orig <- glucose_pick_rows(transports$orig, split$test_days)
    error_table <- glucose_holdout_error_table(pred_orig, obs_orig, loaded$subject_ids, split$test_days, support_grid)
    selected <- utils::head(error_table, top_n)

    sign_class <- vapply(selected$observed_signs, glucose_classify_sign_string, character(1))
    has_opposite <- length(unique(sign_class)) == 2 && all(c("positive", "negative") %in% sign_class)
    all_match <- all(selected$observed_signs == selected$predicted_signs)

    signed_means <- numeric(nrow(selected))
    for (k in seq_len(nrow(selected))) {
      subject_index <- selected$subject_index[k]
      days_k <- as.numeric(strsplit(selected$holdout_days[k], ",", fixed = TRUE)[[1]])
      obs_mat <- obs_orig[[subject_index]][as.character(days_k), , drop = FALSE]
      disp_mat <- obs_mat - matrix(rep(support_grid, each = nrow(obs_mat)), nrow = nrow(obs_mat))
      signed_means[k] <- mean(trapz_rows(support_grid, disp_mat))
    }

    contrast <- if (length(signed_means) >= 2) abs(diff(sort(signed_means)))[1] else 0

    rows[[length(rows) + 1]] <- data.frame(
      seed = seed,
      top1_index = selected$subject_index[1],
      top1_id = selected$subject_id[1],
      top1_days = selected$holdout_days[1],
      top1_err = selected$mean_holdout_w1[1],
      top1_sign = sign_class[1],
      top2_index = selected$subject_index[2],
      top2_id = selected$subject_id[2],
      top2_days = selected$holdout_days[2],
      top2_err = selected$mean_holdout_w1[2],
      top2_sign = sign_class[2],
      total_err = sum(selected$mean_holdout_w1),
      contrast = contrast,
      has_opposite = has_opposite,
      all_match = all_match,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  out <- out[order(-out$has_opposite, -out$all_match, out$total_err, -out$contrast), ]
  out$rank <- seq_len(nrow(out))
  out[, c(
    "rank",
    "seed",
    "top1_index",
    "top1_id",
    "top1_days",
    "top1_err",
    "top1_sign",
    "top2_index",
    "top2_id",
    "top2_days",
    "top2_err",
    "top2_sign",
    "total_err",
    "contrast",
    "has_opposite",
    "all_match"
  )]
}
