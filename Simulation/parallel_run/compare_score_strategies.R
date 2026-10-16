#!/usr/bin/env Rscript
# Standalone previews; this script does not modify manuscript files.
# Run from any directory: Rscript /path/to/compare_score_strategies.R
# The default is the published latent_CE_vs_IN comparison.
# Optional matched audit: matched old_root=PATH results_dir=PATH outdir=PATH.
# Criterion: mean_i integral_t [integral_u |estimated - true transport| du]^2 dt.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_dir <- dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg), fixed = TRUE)))
args <- commandArgs(trailingOnly = TRUE)
mode_arg <- args[!grepl("=", args, fixed = TRUE)]
if (length(mode_arg) > 1L || (length(mode_arg) && !mode_arg %in% c("latent_CE_vs_IN", "matched"))) {
  stop("Optional mode: latent_CE_vs_IN (default) or matched")
}
kv_args <- args[grepl("=", args, fixed = TRUE)]
kv <- setNames(as.list(sub("^[^=]*=", "", kv_args)), sub("=.*$", "", kv_args))
if (length(setdiff(names(kv), c("results_dir", "old_root", "outdir")))) stop("Unknown argument.")
mixed_scores <- !length(mode_arg) || mode_arg == "latent_CE_vs_IN"
if (!mixed_scores && is.null(kv$old_root)) {
  stop("The optional matched audit requires old_root=PATH containing setting1/2/{proposed,baselines}; these alternative runs are not required by the paper workflow.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
library(ggplot2)
old_root <- kv$old_root
new_root <- if (is.null(kv$results_dir)) file.path(script_dir, "results", "paper_comparison") else kv$results_dir
outdir <- if (is.null(kv$outdir)) file.path(new_root, "comparison_previews", if (mixed_scores) "latent_CE_vs_IN" else "matched") else kv$outdir
keys <- c("rep", "n", "N", "m")
expected <- expand.grid(rep = 1:200, n = c(20, 50, 100, 200, 400),
                        N = c(3, 5, 10, 20), m = 200)
make_key <- function(d) do.call(paste, c(d[keys], sep = "/"))
expected_keys <- make_key(expected)

read_checked <- function(path, metrics) {
  if (!file.exists(path)) stop("Missing input: ", path)
  d <- readRDS(path)
  if (!is.data.frame(d) || anyDuplicated(names(d))) stop("Invalid table: ", path)
  missing <- setdiff(c(keys, "ok", metrics), names(d))
  if (length(missing)) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  for (key in keys) {
    v <- d[[key]]
    if (!is.numeric(v) || any(!is.finite(v)) || any(v != floor(v))) {
      stop("Invalid key ", key, " in ", path)
    }
  }
  actual_keys <- make_key(d)
  if (anyDuplicated(actual_keys) || !setequal(expected_keys, actual_keys) || nrow(d) != 4000L) {
    stop("Input does not contain exactly the planned 4,000 unique cells/replications: ", path)
  }
  if (!is.logical(d$ok) || anyNA(d$ok) || !all(d$ok)) stop("Failed fits in ", path)
  for (metric in metrics) {
    v <- d[[metric]]
    if (!is.numeric(v) || any(!is.finite(v)) || any(v < 0)) {
      stop("Invalid squared IMSE in ", path, ", field ", metric)
    }
  }
  d[match(expected_keys, actual_keys), c(keys, metrics), drop = FALSE]
}

strategies <- c("IN_vs_IN", "latent_CE_vs_CE")
if (mixed_scores) strategies <- "latent_CE_vs_IN"
strategy_titles <- c(IN_vs_IN = "Proposed IN vs IN benchmarks",
                     latent_CE_vs_CE = "Proposed latent CE vs CE benchmarks",
                     latent_CE_vs_IN = "Proposed latent CE vs IN benchmarks")
all_long <- all_merged <- list()
input_paths <- character()
for (strategy in strategies) {
  for (setting in 1:2) {
    setting_dir <- paste0("setting", setting)
    prop_dir <- if (strategy == "IN_vs_IN") file.path(old_root, setting_dir, "proposed") else
      file.path(new_root, "proposed_latent_CE", setting_dir)
    base_dir <- if (strategy %in% c("IN_vs_IN", "latent_CE_vs_IN")) file.path(new_root, "benchmarks_IN", setting_dir) else
      file.path(old_root, setting_dir, "baselines")
    paths <- file.path(c(prop_dir, base_dir), "raw_results.rds")
    input_paths <- c(input_paths, paths)
    prop <- read_checked(paths[1], "imse_sq")
    base <- read_checked(paths[2], c("imse_sq_lqd", "imse_sq_tangent"))
    d <- prop[keys]
    d$strategy <- strategy
    d$setting <- setting
    d$imse_sq_prop <- prop$imse_sq
    d$imse_sq_lqd <- base$imse_sq_lqd
    d$imse_sq_tangent <- base$imse_sq_tangent
    all_merged[[length(all_merged) + 1L]] <- d
    for (method in c("LQD", "Tangent", "Proposed")) {
      column <- c(LQD = "imse_sq_lqd", Tangent = "imse_sq_tangent", Proposed = "imse_sq_prop")[[method]]
      long <- d[c("strategy", "setting", keys)]
      long$method <- method
      long$imse_sq <- d[[column]]
      all_long[[length(all_long) + 1L]] <- long
    }
  }
}
merged <- do.call(rbind, all_merged)
long <- do.call(rbind, all_long)
rownames(merged) <- rownames(long) <- NULL
stopifnot(nrow(merged) == 8000L * length(strategies),
          nrow(long) == 24000L * length(strategies))

summary_rows <- improvement_rows <- paired_rows <- list()
for (strategy in strategies) for (setting in 1:2) {
  for (n0 in c(20, 50, 100, 200, 400)) for (N0 in c(3, 5, 10, 20)) {
    d <- merged[merged$strategy == strategy & merged$setting == setting &
                  merged$n == n0 & merged$N == N0, ]
    stopifnot(nrow(d) == 200L)
    means <- c(Proposed = mean(d$imse_sq_prop), LQD = mean(d$imse_sq_lqd),
               Tangent = mean(d$imse_sq_tangent))
    if (any(means[c("LQD", "Tangent")] <= 0)) stop("Nonpositive benchmark mean")
    for (method in names(means)) {
      column <- c(Proposed = "imse_sq_prop", LQD = "imse_sq_lqd", Tangent = "imse_sq_tangent")[[method]]
      v <- d[[column]]
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        strategy, setting, n = n0, N = N0, m = 200L, method, reps = length(v),
        mean_imse_sq = mean(v), mcse_imse_sq = sd(v) / sqrt(length(v)),
        median_imse_sq = median(v), sd_imse_sq = sd(v), iqr_imse_sq = IQR(v))
    }
    improvement_rows[[length(improvement_rows) + 1L]] <- data.frame(
      strategy, setting, n = n0, N = N0, m = 200L, reps = nrow(d),
      mean_imse_sq_prop = unname(means["Proposed"]),
      mean_imse_sq_lqd = unname(means["LQD"]),
      mean_imse_sq_tangent = unname(means["Tangent"]),
      improvement_lqd_pct = unname(100 * (1 - means["Proposed"] / means["LQD"])),
      improvement_tangent_pct = unname(100 * (1 - means["Proposed"] / means["Tangent"])))
    for (benchmark in c("LQD", "Tangent")) {
      benchmark_column <- c(LQD = "imse_sq_lqd", Tangent = "imse_sq_tangent")[[benchmark]]
      delta <- d$imse_sq_prop - d[[benchmark_column]]
      paired_rows[[length(paired_rows) + 1L]] <- data.frame(
        strategy, setting, n = n0, N = N0, m = 200L, benchmark, reps = length(delta),
        mean_difference_prop_minus_benchmark = mean(delta),
        mcse_paired_difference = sd(delta) / sqrt(length(delta)),
        fraction_reps_proposed_lower = mean(delta < 0))
    }
  }
}
summary <- do.call(rbind, summary_rows)
improvements <- do.call(rbind, improvement_rows)
paired <- do.call(rbind, paired_rows)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary, file.path(outdir, "summary_by_cell.csv"), row.names = FALSE)
write.csv(improvements, file.path(outdir, "table1_improvement.csv"), row.names = FALSE)
write.csv(paired, file.path(outdir, "paired_differences.csv"), row.names = FALSE)
write.csv(merged, file.path(outdir, "matched_squared_results.csv"), row.names = FALSE)
write.csv(long, file.path(outdir, "all_cells_long.csv"), row.names = FALSE)
write.csv(data.frame(path = normalizePath(input_paths), md5 = unname(tools::md5sum(input_paths))),
          file.path(outdir, "input_checksums.csv"), row.names = FALSE)

plot_data <- long[long$N %in% c(3, 5, 10), ]
write.csv(plot_data, file.path(outdir, "figure3_plot_data.csv"), row.names = FALSE)
# Identical limits for the two score strategies within each setting. No outliers are removed.
ymax <- vapply(1:2, function(s) max(plot_data$imse_sq[plot_data$setting == s]) * 1.04, numeric(1))
write.csv(data.frame(setting = 1:2, ymin = 0, ymax), file.path(outdir, "figure_y_ranges.csv"), row.names = FALSE)
fill_cols <- c(LQD = "#d9d9d9", Tangent = "#92c5de", Proposed = "#b2182b")
line_cols <- c(LQD = "#8c8c8c", Tangent = "#2166ac", Proposed = "#67001f")
for (strategy in strategies) {
  strategy_dir <- if (mixed_scores) outdir else file.path(outdir, strategy)
  dir.create(strategy_dir, showWarnings = FALSE)
  tab <- improvements[improvements$strategy == strategy, ]
  write.csv(tab, file.path(strategy_dir, "table1_improvement.csv"), row.names = FALSE)
  table_rows <- c("% Percentage improvement: 100*(1 - proposed mean / benchmark mean).",
                  "% Negative numbers mean higher error for the proposed method.",
                  "% Columns: Setting, n, LQD N=3,5,10,20, Tangent N=3,5,10,20.")
  for (setting in 1:2) {
    if (setting == 2L) table_rows <- c(table_rows, "\\midrule")
    table_rows <- c(table_rows, sprintf("\\multirow{5}{*}{Setting %d}", setting))
    for (n0 in c(20, 50, 100, 200, 400)) {
      d <- tab[tab$setting == setting & tab$n == n0, ]
      d <- d[match(c(3, 5, 10, 20), d$N), ]
      table_rows <- c(table_rows, paste0("& ", n0, " & ",
        paste(sprintf("%.2f", c(d$improvement_lqd_pct, d$improvement_tangent_pct)), collapse = " & "), " \\\\"))
    }
    pdat <- plot_data[plot_data$strategy == strategy & plot_data$setting == setting, ]
    stopifnot(nrow(pdat) == 9000L)
    write.csv(pdat, file.path(strategy_dir, sprintf("figure3_setting%d_plot_data.csv", setting)), row.names = FALSE)
    pdat$n <- factor(pdat$n, levels = c(20, 50, 100, 200, 400))
    pdat$N <- factor(pdat$N, levels = c(3, 5, 10), labels = paste0("N = ", c(3, 5, 10)))
    pdat$method <- factor(pdat$method, levels = c("LQD", "Tangent", "Proposed"))
    legend_labels <- if (mixed_scores) c(LQD = "LQD (IN)", Tangent = "Tangent (IN)",
                                        Proposed = "Proposed (latent CE)") else
      c(LQD = "LQD", Tangent = "Tangent", Proposed = "Proposed")
    p <- ggplot(pdat, aes(x = n, y = imse_sq)) +
      geom_boxplot(aes(fill = method, colour = method),
        position = position_dodge2(width = .95, preserve = "single", padding = .05),
        width = .85, outlier.size = .7, linewidth = .45) +
      facet_wrap(~ N, nrow = 1, drop = FALSE) +
      scale_fill_manual(values = fill_cols, labels = legend_labels) +
      scale_colour_manual(values = line_cols, guide = "none") +
      scale_y_continuous(limits = c(0, ymax[setting]), expand = expansion(mult = 0)) +
      labs(x = "Sample size n", y = "IMSE", fill = "Method",
           title = sprintf("Setting %d: %s", setting, strategy_titles[[strategy]])) +
      theme_bw(base_size = 12) +
      theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 14),
        legend.position = "bottom", legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14), legend.key.size = grid::unit(1.1, "cm"),
        legend.spacing.x = grid::unit(.4, "cm"),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(face = "bold", size = 14), panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = .5, face = "bold", size = 15)) +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE))
    for (extension in c("pdf", "png")) {
      path <- file.path(strategy_dir, sprintf("Figure3_setting%d.%s", setting, extension))
      ggsave(path, p, width = 12, height = 4.8, dpi = 300)
    }
  }
  writeLines(table_rows, file.path(strategy_dir, "table1_rows.tex"))
}
writeLines(c(
  sprintf("PASS: all %d input tables contain exactly the planned 4,000 unique (rep,n,N,m) keys.", length(input_paths)),
  "All planned rows have ok=TRUE and finite nonnegative squared IMSE; no rows were dropped.",
  "Each strategy has 40 cells: 2 settings x 5 sample sizes x 4 observation frequencies, with 200 repetitions each.",
  if (mixed_scores) "latent_CE_vs_IN: proposed latent=TRUE, methodXi=CE; both benchmark methods use IN." else c(
    "IN_vs_IN: proposed latent=FALSE, methodXi=IN; both benchmark methods use IN.",
    "latent_CE_vs_CE: proposed latent=TRUE, methodXi=CE; both benchmark methods use CE."),
  "Root validation separately verifies estimator metadata and agreement with saved squared-Wasserstein curves.",
  "Table 1: 100*(1 - mean proposed IMSE / mean benchmark IMSE); positive is improvement, negative is deterioration.",
  "Figure 3 previews show N=3,5,10 and all 200 repetitions per cell; all outliers are displayed.",
  if (mixed_scores) "Each setting's vertical limits include the largest observed value; settings may have different limits." else
    "Vertical limits are common between strategies within each setting and include the largest observed value.",
  "No manuscript, supplement, response or existing results were edited."), file.path(outdir, "validation_report.txt"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("Wrote comparison previews to", normalizePath(outdir), "\n")
