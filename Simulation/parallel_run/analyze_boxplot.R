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
# ============================================================

rm(list = ls())

# -----------------------------
# USER PARAMETERS: edit by hand
# -----------------------------

# Option 1: use merged_raw_results.rds produced by your analysis script
# merged_path <- "plots/LatentF_CE_095_8_b1_beta063/merged_raw_results.rds"

# Option 2: if you do NOT want to use merged_path, set merged_path <- NULL
# and provide prop_dir + base_dir below
# prop_dir <- NULL
# base_dir <- NULL
prop_dir <- "Res/Res_LatentF_CE_095_8_b1_beta063"
base_dir <- "Res/Res_baselines_mse_t_b1_beta063"

# Output folder
outdir <- "plots/LatentF_CE_095_8_b1_beta063/figures"

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
read_merged_results <- function(merged_path = NULL, prop_dir = NULL, base_dir = NULL) {
  
  # Case 1: read merged file directly
  if (!is.null(merged_path) && file.exists(merged_path)) {
    merged <- readRDS(merged_path)
    
    req <- c("rep", "n", "N", "m", "imse_prop", "imse_lqd", "imse_tangent")
    miss <- setdiff(req, names(merged))
    if (length(miss) > 0) {
      stop("merged_path exists but is missing required columns: ",
           paste(miss, collapse = ", "))
    }
    
    return(merged[, req, drop = FALSE])
  }
  
  # Case 2: merge prop/base raw results
  if (is.null(prop_dir) || is.null(base_dir)) {
    stop("Please either:\n",
         "  1) provide merged_path, OR\n",
         "  2) provide both prop_dir and base_dir.")
  }
  
  prop_raw_path <- file.path(prop_dir, "raw_results.rds")
  base_raw_path <- file.path(base_dir, "raw_results.rds")
  
  if (!file.exists(prop_raw_path)) stop("Missing: ", prop_raw_path)
  if (!file.exists(base_raw_path)) stop("Missing: ", base_raw_path)
  
  prop_raw <- readRDS(prop_raw_path)
  base_raw <- readRDS(base_raw_path)
  
  if (!("ok" %in% names(prop_raw))) stop("prop raw_results.rds must contain column 'ok'")
  if (!("ok" %in% names(base_raw))) stop("base raw_results.rds must contain column 'ok'")
  
  prop_ok <- prop_raw[prop_raw$ok, , drop = FALSE]
  base_ok <- base_raw[base_raw$ok, , drop = FALSE]
  
  required_prop <- c("rep", "n", "N", "m", "imse")
  required_base <- c("rep", "n", "N", "m", "imse_lqd", "imse_tangent")
  
  miss_prop <- setdiff(required_prop, names(prop_ok))
  miss_base <- setdiff(required_base, names(base_ok))
  
  if (length(miss_prop) > 0) {
    stop("prop raw_results.rds missing required columns: ",
         paste(miss_prop, collapse = ", "))
  }
  if (length(miss_base) > 0) {
    stop("base raw_results.rds missing required columns: ",
         paste(miss_base, collapse = ", "))
  }
  
  prop_ok <- prop_ok[, c("rep", "n", "N", "m", "imse"), drop = FALSE]
  names(prop_ok)[names(prop_ok) == "imse"] <- "imse_prop"
  
  base_ok <- base_ok[, c("rep", "n", "N", "m", "imse_lqd", "imse_tangent"), drop = FALSE]
  
  merged <- merge(prop_ok, base_ok, by = c("rep", "n", "N", "m"), all = FALSE)
  
  if (nrow(merged) == 0) {
    stop("After merging proposed and baseline results, there are zero rows.\n",
         "Check whether (rep, n, N, m) match between the two runs.")
  }
  
  merged
}

make_long_df <- function(merged_df) {
  lqd_df <- data.frame(
    rep = merged_df$rep,
    n = merged_df$n,
    N = merged_df$N,
    m = merged_df$m,
    method = "LQD",
    imse = merged_df$imse_lqd
  )
  
  tan_df <- data.frame(
    rep = merged_df$rep,
    n = merged_df$n,
    N = merged_df$N,
    m = merged_df$m,
    method = "Tangent",
    imse = merged_df$imse_tangent
  )
  
  prop_df <- data.frame(
    rep = merged_df$rep,
    n = merged_df$n,
    N = merged_df$N,
    m = merged_df$m,
    method = "Proposed",
    imse = merged_df$imse_prop
  )
  
  out <- rbind(lqd_df, tan_df, prop_df)
  out <- out[is.finite(out$imse), , drop = FALSE]
  out
}

# -----------------------------
# READ DATA
# -----------------------------
merged <- read_merged_results(
  merged_path = merged_path,
  prop_dir = prop_dir,
  base_dir = base_dir
)

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
  
  p <- ggplot(long_df, aes(x = n, y = imse)) +
    box_layer +
    facet_wrap(~ N, nrow = 1, drop = FALSE) +
    scale_fill_manual(values = fill_cols) +
    scale_color_manual(values = line_cols, guide = "none") +
    labs(
      x = "Sample size n",
      y = "IMSE",
      fill = "Method"
    ) + ylim(c(0.02,0.1))+
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
  
  print(p)
  
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