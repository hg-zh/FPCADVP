#!/usr/bin/env Rscript
## plot_phase_transition.R
## ------------------------------------------------------------
## Figure 2, left: direct-U FPCA with numerical-integration scores.
## Optional CLI arguments: outdir=/path/to/results, display=paper (default) or full.
## ------------------------------------------------------------
options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
source_file <- if (length(sys.frames())) sys.frames()[[1]]$ofile else NULL
script_path <- if (length(script_arg)) gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE) else source_file
if (is.null(script_path)) {
  candidates <- c("plot_phase_transition.R", "Code/Simulation/Rates_Figure2/plot_phase_transition.R")
  script_path <- candidates[file.exists(candidates)][1]
}
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))

suppressPackageStartupMessages(library(ggplot2))

read_summary <- function(path) {
  if (!file.exists(path)) stop("Cannot find summary file: ", path)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

display_curve_slopes <- function(df) {
  out <- lapply(sort(unique(df$N)), function(N_val) {
    dd <- df[df$N == N_val, , drop = FALSE]
    
    fit <- stats::lm(log(mean_imse) ~ log(n), data = dd)
    
    data.frame(
      N = N_val,
      n_min = min(dd$n),
      n_max = max(dd$n),
      displayed_slope = unname(stats::coef(fit)[["log(n)"]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

## ------------------------------------------------------------
## User settings
## ------------------------------------------------------------
outdir <- "outputs/estimator13_IN"
outdir_arg <- grep("^outdir=", commandArgs(trailingOnly = TRUE), value = TRUE)
if (length(outdir_arg)) outdir <- sub("^outdir=", "", tail(outdir_arg, 1))
if (!grepl("^/", outdir)) outdir <- file.path(script_dir, outdir)
summary_path <- file.path(outdir, "summary_densegrid.csv")
display_arg <- grep("^display=", commandArgs(trailingOnly = TRUE), value = TRUE)
display <- if (length(display_arg)) sub("^display=", "", tail(display_arg, 1)) else "paper"
if (!display %in% c("paper", "full")) stop("display must be 'paper' or 'full'.")
file_stem <- if (display == "paper") "rate_decay_IN" else "rate_decay_IN_full_grid"
outfile <- file.path(outdir, "figures", paste0(file_stem, ".pdf"))

m_target <- 200
slope_n_min <- 70
show_error_bars <- FALSE

## theory line placement
anchor_N <- 10
theory_shift <- -0.04

## ------------------------------------------------------------
## Read and prepare data
## ------------------------------------------------------------
sum_df <- read_summary(summary_path)

need_cols <- c("n", "N", "m", "mean_imse", "a_poly", "theoretical_slope")
missing_cols <- setdiff(need_cols, names(sum_df))
if (length(missing_cols) > 0) {
  stop("summary file is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

df <- sum_df[sum_df$m == m_target, , drop = FALSE]
if (nrow(df) == 0) stop("No rows found for m = ", m_target)
if (any(!is.finite(df$mean_imse) | df$mean_imse <= 0 |
        !is.finite(df$n) | df$n <= 0 | !is.finite(df$N) | df$N <= 0)) {
  stop("All summary cells must have finite positive n, N, and mean_imse.")
}
if (anyDuplicated(df[c("n", "N")])) stop("Duplicate (n, N) cells in summary.")
if (nrow(df) != length(unique(df$n)) * length(unique(df$N))) {
  stop("Summary does not contain the complete (n, N) grid.")
}

# The manuscript display uses the requested range. Preserve the full simulation
# summary and provide display=full to reproduce all 42 simulated cells.
total_cells <- nrow(df)
if (display == "paper") {
  n_display <- c(20L, 30L, 40L, 50L, 70L, 100L)
  N_display <- c(3L, 4L, 5L, 7L, 10L)
  df <- df[df$n %in% n_display & df$N %in% N_display, , drop = FALSE]
  if (nrow(df) != length(n_display) * length(N_display)) {
    stop("Summary does not contain the complete manuscript display grid.")
  }
}

df$log_n <- log(df$n)
df$log_imse <- log(df$mean_imse)
df <- df[order(df$N, df$n), , drop = FALSE]
N_levels <- sort(unique(df$N))

if (show_error_bars && "se_imse" %in% names(df)) {
  df$imse_lo <- pmax(df$mean_imse - 1.96 * df$se_imse, .Machine$double.eps)
  df$imse_hi <- df$mean_imse + 1.96 * df$se_imse
  df$log_lo <- log(df$imse_lo)
  df$log_hi <- log(df$imse_hi)
}

## Descriptive slopes of log(mean IMSE) over the exact displayed range.
slope_df <- display_curve_slopes(df)
print(slope_df)

## ------------------------------------------------------------
## Legend labels: just N values
## ------------------------------------------------------------
line_labels <- as.character(N_levels)
df$series <- factor(as.character(df$N), levels = line_labels)

## ------------------------------------------------------------
## Theoretical slope line
## ------------------------------------------------------------
a_poly <- unique(df$a_poly)
if (length(a_poly) != 1 || !is.finite(a_poly) || a_poly <= 1) {
  stop("Summary must specify one finite eigenvalue-decay exponent a > 1.")
}
# Corollary 1's algebraic benchmark, omitting logarithmic factors.
# This is -1/8 = -0.125 when a = 2; the vertical placement is arbitrary.
ref_slope <- -(a_poly - 1) / (3 * a_poly + 2)
if (any(!is.finite(df$theoretical_slope)) ||
    any(abs(df$theoretical_slope - ref_slope) > 1e-10)) {
  stop("Saved reference slope disagrees with -(a-1)/(3a+2); rerun the corrected driver.")
}

anchor_df <- df[df$N == anchor_N & df$n >= slope_n_min, , drop = FALSE]
if (nrow(anchor_df) == 0) anchor_df <- df[df$N == anchor_N, , drop = FALSE]
if (nrow(anchor_df) == 0) stop("No rows found for anchor_N = ", anchor_N)
x0 <- mean(anchor_df$log_n)
y0 <- mean(anchor_df$log_imse) + theory_shift
x_pad <- max(0.025 * diff(range(df$log_n)), 0.025)
xlim <- range(df$log_n) + c(-x_pad, x_pad)
x_ref <- seq(xlim[1], xlim[2], length.out = 200)
ref_df <- data.frame(log_n = x_ref, log_imse = y0 + ref_slope * (x_ref - x0))
message(sprintf("Algebraic reference slope: %.6f; displaying %d of %d cells, n = %s, N = %s.",
                ref_slope, nrow(df), total_cells, paste(sort(unique(df$n)), collapse = ", "),
                paste(N_levels, collapse = ", ")))

## ------------------------------------------------------------
## Show every cell in the documented display range.
## ------------------------------------------------------------
df_plot <- df

if (show_error_bars && all(c("log_lo", "log_hi") %in% names(df_plot))) {
  ylim <- range(c(df_plot$log_lo, df_plot$log_hi, df_plot$log_imse), finite = TRUE)
} else {
  ylim <- range(df_plot$log_imse, finite = TRUE)
}

if (!is.null(ref_df) && nrow(ref_df) > 0) {
  ylim <- range(c(ylim, ref_df$log_imse), finite = TRUE)
}

y_pad <- 0.08 * diff(ylim)
if (!is.finite(y_pad) || y_pad <= 0) y_pad <- 0.05
ylim <- c(ylim[1] - y_pad, ylim[2] + y_pad)

## ------------------------------------------------------------
## Aesthetics
## ------------------------------------------------------------
cols <- grDevices::colorRampPalette(
  c("#D55E00", "#E69F00", "#56B4E9", "#009E73", "#0072B2")
)(length(N_levels))
names(cols) <- line_labels

legend_breaks <- line_labels

## ------------------------------------------------------------
## Build ggplot
## ------------------------------------------------------------
p <- ggplot() +
  geom_line(
    data = df_plot,
    aes(x = log_n, y = log_imse, color = series, group = series),
    linewidth = 1.2
  ) +
  geom_point(
    data = df_plot,
    aes(x = log_n, y = log_imse, color = series),
    size = 2.2,
    show.legend = FALSE
  )

if (show_error_bars && all(c("log_lo", "log_hi") %in% names(df_plot))) {
  p <- p + geom_errorbar(
    data = df_plot,
    aes(x = log_n, y = log_imse, ymin = log_lo, ymax = log_hi, color = series),
    width = 0.02,
    show.legend = FALSE
  )
}

if (!is.null(ref_df)) {
  p <- p + geom_line(
    data = ref_df,
    aes(x = log_n, y = log_imse),
    color = "black",
    linetype = "dashed",
    linewidth = 1.2,
    inherit.aes = FALSE,
    show.legend = FALSE
  )
}

p <- p +
  scale_color_manual(
    values = cols,
    breaks = legend_breaks,
    name = "N"
  ) +
  scale_x_continuous(n.breaks = 5) +
  scale_y_continuous(n.breaks = 4) +
  coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(
    x = "log(n)",
    y = "log(IMSE)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey92"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20),
    legend.position = "right",
    legend.justification = "center",
    legend.margin = margin(l = 6),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.key.width = grid::unit(1, "cm")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 0.9,
        shape = NA
      )
    )
  )

## ------------------------------------------------------------
## Show in RStudio Plots pane
## ------------------------------------------------------------
if (interactive()) print(p)

## ------------------------------------------------------------
## Save without opening an incidental Rplots.pdf in noninteractive runs.
## ------------------------------------------------------------
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
ggsave(outfile, plot = p, width = 500/72, height = 400/72, units = "in", bg = "white")
write.csv(df[c("n", "N", "m", "reps_ok", "mean_imse", "se_imse", "theoretical_slope")],
          file.path(dirname(outfile), paste0(file_stem, "_cells.csv")), row.names = FALSE)
write.csv(slope_df, file.path(dirname(outfile), paste0(file_stem, "_slopes.csv")), row.names = FALSE)
writeLines(c(
  paste("Display:", display),
  sprintf("Shown cells: %d of %d; n = %s; N = %s.", nrow(df), total_cells,
          paste(sort(unique(df$n)), collapse = ", "), paste(N_levels, collapse = ", ")),
  "The manuscript display excludes N=20 and n=150 (log(n) approximately 5).",
  "Complete raw results, summaries and full-range slope estimates are retained unchanged.",
  "Use display=full to reproduce the complete simulation grid.",
  sprintf("Reference slope: %.6f; reference height is arbitrary.", ref_slope)
), file.path(dirname(outfile), paste0(file_stem, "_display.txt")))

