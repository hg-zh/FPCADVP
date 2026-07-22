## plot_phase_transition_ggplot.R
## ------------------------------------------------------------
## Interactive ggplot2 script for RStudio
## ------------------------------------------------------------
setwd("~/Documents/research/Hans/random object/JRSSB_RR1/Final_Code/Simulation/Rates_Figure2")
options(stringsAsFactors = FALSE)

## Close old graphics devices if needed
try(graphics.off(), silent = TRUE)

library(ggplot2)

read_summary <- function(path) {
  if (!file.exists(path)) stop("Cannot find summary file: ", path)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

tail_curve_slopes <- function(df, slope_n_min) {
  out <- lapply(sort(unique(df$N)), function(N_val) {
    dd <- df[df$N == N_val & df$n >= slope_n_min, , drop = FALSE]
    if (nrow(dd) < 2) dd <- df[df$N == N_val, , drop = FALSE]
    
    fit <- stats::lm(log(mean_imse) ~ log(n), data = dd)
    
    data.frame(
      N = N_val,
      tail_slope = unname(stats::coef(fit)[["log(n)"]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

## ------------------------------------------------------------
## User settings
## ------------------------------------------------------------
outdir <- "outputs"
summary_path <- file.path(outdir, "summary_densegrid.csv")
outfile <- file.path(outdir, "figures", "loglog_phase_transition_densegrid.pdf")

m_target <- 200
slope_n_min <- 70
show_error_bars <- FALSE

## visible x-range
x_min_plot <- 2.95
x_max_plot <- 5.05

## theory line placement
anchor_N <- 10
theory_shift <- -0.04

## ------------------------------------------------------------
## Read and prepare data
## ------------------------------------------------------------
sum_df <- read_summary(summary_path)

need_cols <- c("n", "N", "m", "mean_imse")
missing_cols <- setdiff(need_cols, names(sum_df))
if (length(missing_cols) > 0) {
  stop("summary file is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

df <- sum_df[
  sum_df$m == m_target &
    is.finite(sum_df$mean_imse) &
    sum_df$mean_imse > 0,
  ,
  drop = FALSE
]

if (nrow(df) == 0) stop("No finite positive rows found for m = ", m_target)

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

## optional: compute tail slopes for your own checking
slope_df <- tail_curve_slopes(df, slope_n_min = slope_n_min)
print(slope_df)

## ------------------------------------------------------------
## Legend labels: just N values
## ------------------------------------------------------------
line_labels <- as.character(N_levels)
df$series <- factor(as.character(df$N), levels = line_labels)

## ------------------------------------------------------------
## Theoretical slope line
## ------------------------------------------------------------
ref_slope <- if ("theoretical_slope" %in% names(df)) df$theoretical_slope[1] else NA_real_
a_poly <- if ("a_poly" %in% names(df)) df$a_poly[1] else NA_real_

if (is.finite(ref_slope)) {
  anchor_df <- df[df$N == anchor_N & df$n >= slope_n_min, , drop = FALSE]
  if (nrow(anchor_df) == 0) anchor_df <- df[df$N == anchor_N, , drop = FALSE]
  if (nrow(anchor_df) == 0) stop("No rows found for anchor_N = ", anchor_N)
  
  x0 <- mean(anchor_df$log_n)
  y0 <- mean(anchor_df$log_imse) + theory_shift
  
  x_ref <- seq(x_min_plot, x_max_plot, length.out = 200)
  y_ref <- y0 + ref_slope * (x_ref - x0)
  
  ref_df <- data.frame(
    log_n = x_ref,
    log_imse = y_ref
  )
} else {
  ref_df <- NULL
}

## ------------------------------------------------------------
## Keep only visible data
## ------------------------------------------------------------
df_plot <- df[df$log_n >= x_min_plot & df$log_n <= x_max_plot, , drop = FALSE]
if (nrow(df_plot) == 0) {
  stop("No data remain inside the chosen x-range.")
}

xlim <- c(x_min_plot, x_max_plot)

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
  scale_y_continuous(breaks = c(-6.0, -5.5)) +
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
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
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
print(p)

## ------------------------------------------------------------
## Optional: save after you are happy with the plot
## ------------------------------------------------------------
# dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
ggsave(outfile, plot = p, width = 500/72, height = 400/72, units = "in", bg = "white")

