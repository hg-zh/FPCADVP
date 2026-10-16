# fig2_plots_modular.R
# ------------------------------------------------------------
# Reproduces the MIDDLE and RIGHT panels of Figure 2 using the
# new modular functions in transport_fpca_modular.R.
#
# This replaces sim_sparse_illu.R + utlis.R.
#
# Outputs (under this script's outputs directory by default):
#   - Fig2_middle_modular.pdf
#   - Fig2_right_modular.html with its companion dependency directory
#
# Required packages: fdapace, ggplot2, plotly, htmlwidgets
# ------------------------------------------------------------

library(fdapace)
library(ggplot2)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript.")
script_dir <- dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg), fixed = TRUE), mustWork = TRUE))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L || (length(args) && !grepl("^outdir=", args))) stop("Optional argument: outdir=PATH")
outdir <- if (length(args)) sub("^outdir=", "", args) else file.path(script_dir, "outputs")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
source(file.path(script_dir, "transport_fpca_modular.R"))

# -------------------------
# Global settings
# -------------------------
K <- 50
m <- 200
xgrid <- seq(0, 1, length.out = 128)

# Link function used in the manuscript
# g: R -> (-1,1), g_inv: (-1,1) -> R
# (Keep as-is to match the simulation section in the paper.)
g <- function(z) 2 * atan(z) / pi
g_inv <- function(u) tan(pi * u / 2)

# =========================
# 1) MIDDLE panel of Fig 2
# =========================
# "Sparsely observed transports and corresponding estimates" for N=4.
# We follow the old sim_sparse_illu.R behavior: use the TRUE transports
# as the "observed" curves (solid), and plot fitted curves (dashed).

set.seed(2001)  # matches old script: r=1 -> seed=r+2000
n_mid <- 100
N_mid <- 4

sim_mid <- simulate_transport_process(
  n = n_mid,
  N = N_mid,
  m = m,
  xgrid = xgrid,
  K = K,
  t_design = "random",
  g = g,
  seed = 2001
)

# Use TRUE transports as observed (solid lines), like the previous illustration
T_obs_mid <- sim_mid$T_true

fit_mid <- fit_transport_fpca(
  T_list = T_obs_mid,
  t_list = sim_mid$t_list,
  xgrid = xgrid,
  kappa = sim_mid$kappa_true,
  g = g,
  g_inv = g_inv,
  methodXi = "CE",        
  fpca_optns = list(methodMuCovEst = "smooth"),
  verbose = TRUE
)

# Match the old illustration subject index if desired
i_show <- 9
obs_times <- sim_mid$t_list[[i_show]]

# Nearest workGrid indices for each observed time
idx_grid <- sapply(obs_times, function(t) which.min(abs(fit_mid$workGrid - t)))

# Build long-format data for ggplot
labs_time <- sprintf("t=%.2f", obs_times)
cols_time <- c("red", "orange", "purple", "blue")
if (length(cols_time) < length(obs_times)) {
  # fallback if N_mid changes
  cols_time <- rep(cols_time, length.out = length(obs_times))
}
names(cols_time) <- labs_time

DF_mid <- do.call(
  rbind,
  lapply(seq_along(obs_times), function(j) {
    data.frame(
      x = xgrid,
      # subtract identity map; flip sign to match the paper's orientation
      y = -(T_obs_mid[[i_show]][j, ] - xgrid),
      y_hat = -(fit_mid$T_hat[[i_show]][idx_grid[j], ] - xgrid),
      time = labs_time[j]
    )
  })
)
DF_mid$time <- factor(DF_mid$time, levels = labs_time)

p_mid <- ggplot(DF_mid, aes(x = x, y = y, color = time)) +
  geom_line(linewidth = 1.2) +
  geom_line(aes(y = y_hat), linewidth = 1.2, linetype = "dashed") +
  scale_color_manual(values = cols_time) +
  labs(x = "u", y = "T(u)-u", color = NULL) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    
    # legend inside, top-right
    legend.position = "inside",
    legend.position.inside = c(0.98, 1),
    legend.justification = c(1, 1),
    
    # optional: make legend easier to read
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_blank()
  )

p_mid

ggsave(file.path(outdir, "Fig2_middle_modular.pdf"), p_mid, width = 500/72, height = 400/72)
# Save
# try(ggsave("Fig2_middle_modular.png", p_mid, width = 4.6, height = 3.6, dpi = 300), silent = TRUE)
# try(ggsave("Fig2_middle_modular.pdf", p_mid, width = 4.6, height = 3.6), silent = TRUE)

# =========================
# 2) RIGHT panel of Fig 2
# =========================
# True (gray) and estimated (blue) transport surfaces over time.
# Red curves are the transports at the observed times.
# Identity is subtracted; sign flipped to match the middle panel: z = -(T(u)-u).

for (package in c("plotly", "htmlwidgets")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Install required package: ", package)
}
library(plotly)
library(htmlwidgets)

# Use the same subject as the middle panel
i_show <- 9

t_grid <- fit_mid$workGrid
u_grid <- xgrid

# ---- True surface on the same time grid (workGrid) ----
# Requires sim_mid$Xi and sim_mid$T0_true from simulate_transport_process().
T_true_grid <- true_transports_on_grid(
  Xi_row = sim_mid$Xi[i_show, ],
  T0_row = sim_mid$T0_true[i_show, ],
  t_grid = t_grid,
  xgrid  = u_grid,
  K      = K,
  g      = g
)

# ---- Estimated surface already on workGrid ----
T_hat_grid <- fit_mid$T_hat[[i_show]]   # dim: length(t_grid) x length(u_grid)

# ---- Build z matrices for plotly (must be length(u) x length(t)) ----
# identity-subtracted & sign-flipped (same as middle plot):
# z(u,t) = -(T(t,u) - u)
Z_true <- -sweep(t(T_true_grid), 1, u_grid, "-")  # dim: G x L
Z_hat  <- -sweep(t(T_hat_grid),  1, u_grid, "-")  # dim: G x L

fig <- plotly::plot_ly()

# Draw TRUE first (gray), then EST (blue) so blue is visually "on top"
fig <- plotly::add_surface(
  fig,
  x = t_grid, y = u_grid, z = Z_true,
  showscale = FALSE,
  opacity = 0.60,
  colorscale = list(list(0, "gray"), list(1, "gray"))
)

fig <- plotly::add_surface(
  fig,
  x = t_grid, y = u_grid, z = Z_hat,
  showscale = FALSE,
  opacity = 0.60,
  colorscale = list(list(0, "blue"), list(1, "blue"))
)

# ---- Add red curves at observed times ----
t_obs <- sim_mid$t_list[[i_show]]

for (j in seq_along(t_obs)) {
  z_line <- -(T_obs_mid[[i_show]][j, ] - u_grid)  # your "observed" (here true) transports
  fig <- plotly::add_trace(
    fig,
    x = rep(t_obs[j], length(u_grid)),
    y = u_grid,
    z = z_line,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "red", width = 6),
    showlegend = FALSE
  )
}

scene <- list(
  xaxis = list(title = "time", titlefont = list(size = 18)),
  yaxis = list(title = "u",    titlefont = list(size = 18)),
  zaxis = list(title = "T(u)-u", titlefont = list(size = 18)),
  camera = list(center = list(x = 0, y = 0, z = -0.2),
                eye    = list(x = 1.25, y = 1.25, z = 1.25))
)

fig <- plotly::layout(fig, scene = scene, showlegend = FALSE)

# Save interactive HTML (selfcontained=FALSE avoids needing pandoc)
htmlwidgets::saveWidget(fig, file.path(outdir, "Fig2_right_modular.html"), selfcontained = FALSE)

# Portable static companion from the same fitted/true arrays. The archived
# manuscript image uses the interactive plot's camera/rendering; this PDF is
# provided as a deterministic, browser-independent view of those same values.
pdf(file.path(outdir, "Fig2_right_static.pdf"), width = 7, height = 5)
view <- persp(t_grid, u_grid, t(Z_hat), theta = 40, phi = 25,
              xlab = "time", ylab = "u", zlab = "T(u)-u",
              zlim = range(c(Z_true, Z_hat)),
              col = grDevices::adjustcolor("blue", alpha.f = 0.35),
              border = NA, ticktype = "detailed")
for (k in unique(round(seq(1, length(u_grid), length.out = 16)))) {
  lines(trans3d(t_grid, rep(u_grid[k], length(t_grid)), Z_true[k, ], view),
        col = grDevices::adjustcolor("grey40", alpha.f = 0.6), lwd = 0.7)
}
for (j in seq_along(t_obs)) {
  z_line <- -(T_obs_mid[[i_show]][j, ] - u_grid)
  lines(trans3d(rep(t_obs[j], length(u_grid)), u_grid, z_line, view),
        col = "red", lwd = 1.5)
}
dev.off()

# Display in RStudio Viewer when sourced interactively.
if (interactive()) fig

