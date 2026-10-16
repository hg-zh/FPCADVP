# ============================================================
# plot_boxplots_by_N_RStudio.R
# ------------------------------------------------------------
# Create a 1x3 boxplot array:
#   Left   = N = 3
#   Middle = N = 5
#   Right  = N = 10
#
# Within each panel:
#   x-axis = sample size n
#   y-axis = IMSE
#   colors = methods
#
# Methods:
#   LQD, Tangent, Proposed
#
# The Proposed method is highlighted for better illustration.
# Approved publication comparison: proposed latent=TRUE/methodXi=CE
# (estimator (14)) versus LQD and Tangent with methodXi=IN.
# ============================================================

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE)))
} else normalizePath(getwd())
args <- commandArgs(trailingOnly = TRUE)
kv <- list()
for (arg in args) {
  pair <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(pair) != 2 || !pair[1] %in% c("prop_dir", "base_dir", "outdir", "setting", "results_dir")) {
    stop("Expected setting=1|2, results_dir=PATH, prop_dir=PATH, base_dir=PATH or outdir=PATH")
  }
  kv[[pair[1]]] <- pair[2]
}
source(file.path(script_dir, "summarize_squared_imse.R"))
result_dir <- if (is.null(kv$results_dir)) file.path(script_dir, "results", "paper_comparison") else kv$results_dir
setting <- if (is.null(kv$setting)) 1L else as.integer(kv$setting)
if (is.na(setting) || !setting %in% c(1L, 2L)) stop("setting must be 1 or 2")

# -----------------------------
# DEFAULT PLOT PARAMETERS
# -----------------------------

# The inputs must contain explicitly named squared-W1 outcomes.
prop_dir <- if (!is.null(kv$prop_dir)) kv$prop_dir else
  file.path(result_dir, "proposed_latent_CE", paste0("setting", setting))
base_dir <- if (!is.null(kv$base_dir)) kv$base_dir else
  file.path(result_dir, "benchmarks_IN", paste0("setting", setting))
outdir <- if (!is.null(kv$outdir)) kv$outdir else file.path(result_dir, "figures", paste0("setting", setting))

# Which m values to plot
# Example: m_values <- 200
#          m_values <- c(50, 200)
#          m_values <- NULL   # plot all available m
m_values <- 200

# Fixed panel order: left to right = N=3, N=5, N=10
N_panels <- c(3, 5, 10)

# Output figure format: "pdf" or "png"
fig_format <- "pdf"

# Figure size
fig_width  <- 12
fig_height <- 4.5
fig_dpi    <- 400

# Whether to use log scale on y-axis
y_log <- FALSE

# Whether to show outliers
show_outliers <- TRUE

# Whether to add a main title
use_main_title <- F

# -----------------------------
# LOAD PACKAGE
# -----------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required but not installed.")
}
library(ggplot2)

# -----------------------------
# HELPER FUNCTIONS
# -----------------------------
make_long_df <- function(merged_df) {
  lqd_df <- data.frame(
    rep = merged_df$rep,
    n = merged_df$n,
    N = merged_df$N,
    m = merged_df$m,
    method = "LQD",
    imse_sq = merged_df$imse_sq_lqd
  )
  
  tan_df <- data.frame(
    rep = merged_df$rep,
    n = merged_df$n,
    N = merged_df$N,
    m = merged_df$m,
    method = "Tangent",
    imse_sq = merged_df$imse_sq_tangent
  )
  
  prop_df <- data.frame(
    rep = merged_df$rep,
    n = merged_df$n,
    N = merged_df$N,
    m = merged_df$m,
    method = "Proposed",
    imse_sq = merged_df$imse_sq_prop
  )
  
  out <- rbind(lqd_df, tan_df, prop_df)
  out
}

# -----------------------------
# READ DATA
# -----------------------------
merged <- read_squared_results(prop_dir, base_dir, setting = setting)

# Filter by m
if (!is.null(m_values)) {
  merged <- merged[merged$m %in% m_values, , drop = FALSE]
}

# Filter by N panels
merged <- merged[merged$N %in% N_panels, , drop = FALSE]

if (nrow(merged) == 0) {
  stop("No rows remain after filtering by m_values / N_panels.")
}

# Create output folder
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Available m values after filtering
m_vals <- sort(unique(merged$m))
cat("Available m values to plot:", paste(m_vals, collapse = ", "), "\n")

# -----------------------------
# COLOR CHOICES
# Proposed is made visually strongest
# -----------------------------
fill_cols <- c(
  "LQD"      = "#d9d9d9",
  "Tangent"  = "#92c5de",
  "Proposed" = "#b2182b"
)

line_cols <- c(
  "LQD"      = "#8c8c8c",
  "Tangent"  = "#2166ac",
  "Proposed" = "#67001f"
)

# -----------------------------
# MAKE FIGURES
# -----------------------------
for (m0 in m_vals) {
  
  d0 <- merged[merged$m == m0, , drop = FALSE]
  if (nrow(d0) == 0) next
  
  long_df <- make_long_df(d0)
  write.csv(long_df, file.path(outdir, sprintf("figure3_squared_plot_data_m%s.csv", m0)),
            row.names = FALSE)
  
  n_levels <- sort(unique(long_df$n))
  long_df$n <- factor(long_df$n, levels = n_levels)
  
  # Fix panel order exactly as N = 3, 5, 10
  long_df$N <- factor(long_df$N, levels = N_panels,
                      labels = paste0("N = ", N_panels))
  
  # Method order inside each n
  long_df$method <- factor(long_df$method, levels = c("LQD", "Tangent", "Proposed"))
  
  if (show_outliers) {
    box_layer <- geom_boxplot(
      aes(fill = method, color = method),
      position = position_dodge2(width = 0.95, preserve = "single", padding = 0.05),
      width = 0.85,
      outlier.size = 0.7,
      linewidth = 0.45
    )
  } else {
    box_layer <- geom_boxplot(
      aes(fill = method, color = method),
      position = position_dodge2(width = 0.95, preserve = "single", padding = 0.05),
      width = 0.85,
      outlier.shape = NA,
      linewidth = 0.45
    )
  }
  
  p <- ggplot(long_df, aes(x = n, y = imse_sq)) +
    box_layer +
    facet_wrap(~ N, nrow = 1, drop = FALSE) +
    scale_fill_manual(values = fill_cols) +
    scale_color_manual(values = line_cols, guide = "none") +
    labs(
      x = "Sample size n",
      y = "IMSE",
      fill = "Method"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 14),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 14),
      legend.text  = element_text(size = 14),
      legend.key.size = grid::unit(1.1, "cm"),
      legend.spacing.x = grid::unit(0.4, "cm"),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold", size = 14),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  
  if (use_main_title) {
    p <- p + ggtitle(sprintf("IMSE boxplots by sample size and method (m = %s)", m0))
  }
  
  if (isTRUE(y_log)) {
    p <- p + scale_y_log10()
  }
  
  if (interactive()) print(p)
  
  outfile <- file.path(
    outdir,
    sprintf("Fig_boxplot_byN_methods_m%s.%s", m0, fig_format)
  )
  
  ggsave(
    filename = outfile,
    plot = p,
    width = fig_width,
    height = fig_height,
    dpi = fig_dpi
  )
  
  cat("Wrote:", normalizePath(outfile), "\n")
}

cat("Done.\n")
