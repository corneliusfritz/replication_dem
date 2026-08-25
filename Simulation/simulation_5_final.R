# Simulation Study 5: Approximation Accuracy vs. Grid Changepoint Steps (Q)
# This script evaluates how the resolution of time-changepoint grid cells (Q) affects
# parameter coverage probability and per-iteration computation times under Blockwise estimation.
# It simulates event sequences, estimates models across a grid of cell sizes, and plots trend curves.

rm(list = ls())
source("Simulation/helper.R")
library(stringr)
library(readr)
library(stringi)
library(redeem)
library(parallel)
library(latex2exp)
library(ggpubr)
library(patchwork)
library(data.table)

# Set random seed to ensure exact reproducibility
set.seed(1234)

# --- 1. Simulation Parameter Setup (N = 500) ---
n_nodes <- 500
continuous_cov <- rnorm(n = n_nodes)
continuous_cov <- outer(X = continuous_cov, Y = continuous_cov, FUN = function(x, y) {
  abs(x - y)
})
categorical_cov <- sample(x = 1:3, size = n_nodes, replace = TRUE)
categorical_cov <- outer(X = categorical_cov, Y = categorical_cov, FUN = function(x, y) {
  x == y
})
popularity_0_1 <- rnorm(n = n_nodes, sd = 1) - (6 + 0.1 * log(n_nodes))
popularity_1_0 <- rnorm(n = n_nodes, sd = 1) - (1.6 + 0.1 * log(n_nodes))

time_changepoints <- seq(0, 10000, length.out = 10)[-c(1, 10)]

baseline_0_1_gt <- -seq(from = 0, to = -0.1, length.out = length(time_changepoints) + 1)[-1]
baseline_1_0_gt <- seq(from = 0, to = -0.1, length.out = length(time_changepoints) + 1)[-1]

# Specify formulas and true coefficients
formula_0_1 <- ~ degrees + current_common_partners(transformation = "log") +
  dyadic_cov(data = continuous_cov) + dyadic_cov(data = categorical_cov)
coef_0_1 <- c(
  "current_common_partners" = -0.5,
  "continuous_cov" = 1,
  "categorical_cov" = 0.5
)
formula_1_0 <- ~ degrees + inertia(transformation = "log") +
  dyadic_cov(data = continuous_cov) + dyadic_cov(data = categorical_cov)
coef_1_0 <- c(
  "inertia" = 0.5,
  "continuous_cov" = 0.5, "categorical_cov" = 0.5
)


simultaneous_interactions <- TRUE
it_max <- 1500
accelerated <- FALSE
estimate_newton <- FALSE
x <- 7
K <- 100
verbose <- FALSE
clust <- makeForkCluster(K, outfile = "outfile.txt")
clusterCall(clust, function() {
  library(redeem)
})
grid_cells <- round(exp(seq(from = 1, to = 7, length.out = 20)) + 1)


res_simulation <- parLapply(
  cl = clust, X = 1:(10 * K), fun = simulation_baseline,
  formula_0_1 = formula_0_1,
  formula_1_0 = formula_1_0,
  coef_0_1 = coef_0_1,
  coef_1_0 = coef_1_0,
  time_changepoints = time_changepoints,
  baseline_0_1_gt = baseline_0_1_gt,
  baseline_1_0_gt = baseline_1_0_gt,
  n_nodes = n_nodes,
  popularity_0_1 = popularity_0_1,
  popularity_1_0 = popularity_1_0,
  simultaneous_interactions = simultaneous_interactions,
  verbose = F,
  it_max = 1000,
  accelerated = FALSE,
  grid_cells = grid_cells
)

saveRDS(res_simulation, file = "Simulation/Results/simulation_5.RDS")
res_simulation <- readRDS("Simulation/Results/simulation_5.RDS")
file.remove("outfile.txt")
stopCluster(clust)

time <- do.call("rbind", lapply(res_simulation, function(x) unlist(x$needed_time)))
colnames(time) <- res_simulation[[1]]$grid_cells
time <- as.data.table(time)

alpha <- 0.05
boot_mean <- function(x) {
  resample <- sample(x, size = length(x), replace = TRUE)
  return(mean(resample))
}
B <- 1000
time_data <- melt.data.table(time)
time_data$value <- time_data$value * 60
time_data <- time_data[, .(
  lower = quantile(value, 0.025),
  min = min(value),
  mean = mean(value),
  ci_lower = quantile(replicate(B, boot_mean(value)), alpha / 2),
  ci_upper = quantile(replicate(B, boot_mean(value)), 1 - alpha / 2),
  upper = quantile(value, 0.975),
  max = max(value, 0.975)
), by = variable]


# quantile(time_data[variable == 1098]$value, 0.025)

my_theme <- theme_pubr() + theme(
  legend.background = element_rect(fill = "transparent", color = NA),
  legend.key = element_rect(fill = "transparent", color = NA),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10),
  plot.title = element_text(hjust = 0.5)
)

plot_tmp <- ggplot(time_data, aes(x = as.numeric(as.character(variable)))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey40", alpha = 0.2) +
  # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "grey20", alpha = 0.3) +
  geom_line(aes(y = mean), color = "black", linewidth = 1) +
  geom_point(aes(y = mean), size = 1.5, color = "black") +
  my_theme +
  theme_pubr(base_size = 17) +
  scale_x_continuous(breaks = pretty(as.numeric(as.character(time_data$variable)))) +
  labs(
    x = "Number of steps of f(Q)",
    y = "Time needed in seconds"
  ) +
  labs(caption = "(a) Per-iteration update time for varying steps of f (Q)") +
  theme(
    plot.caption = element_text(hjust = 0.5, size = 17)
  )



grid_cells <- lapply(res_simulation, function(x) x$grid_cells)


cp_0_1_list <- lapply(res_simulation, function(x) x$cp_0_1_list)
cp_1_0_list <- lapply(res_simulation, function(x) x$cp_1_0_list)


res_0_1 <- mapply(
  function(x, y) {
    cbind(t(do.call("cbind", x)), n_change = res_simulation[[1]]$grid_cells, id = y)
  },
  cp_0_1_list, seq_along(cp_0_1_list),
  SIMPLIFY = FALSE
)
res_0_1 <- do.call("rbind", res_0_1)
res_0_1 <- data.frame(res_0_1)
names(res_0_1)[1:3] <- names(res_simulation[[1]]$coef_0_1)

setDT(res_0_1)
res_1_0 <- mapply(
  function(x, y) {
    cbind(t(do.call("cbind", x)), n_change = res_simulation[[1]]$grid_cells, id = y)
  },
  cp_1_0_list, seq_along(cp_1_0_list),
  SIMPLIFY = FALSE
)
res_1_0 <- do.call("rbind", res_1_0)
res_1_0 <- data.frame(res_1_0)
names(res_1_0)[1:3] <- names(res_simulation[[1]]$coef_1_0)
setDT(res_1_0)

plot_data_1_0 <- res_1_0[, .(
  cp_number_interaction = mean(inertia),
  cp_continuous_cov = mean(continuous_cov),
  cp_categorical_cov = mean(categorical_cov)
), by = n_change]


plot_data_0_1 <- res_0_1[, .(
  cp_current_common_partner = mean(current_common_partners),
  cp_continuous_cov = mean(continuous_cov),
  cp_categorical_cov = mean(categorical_cov)
), by = n_change]

plot_data_0_1 <- melt.data.table(plot_data_0_1, id.vars = "n_change")

master_colors <- c(
  "cp_current_common_partner" = "#0072B2",
  "cp_continuous_cov" = "#E69F00",
  "cp_categorical_cov" = "#D55E00",
  "cp_number_interaction" = "#009E73"
)

master_labels <- c(
  "cp_current_common_partner" = expression(s[i * "," * j * "," * CCP]),
  "cp_continuous_cov" = expression(s[i * "," * j * "," * bold(x)]),
  "cp_categorical_cov" = expression(s[i * "," * j * "," * bold(y)]),
  "cp_number_interaction" = expression(s[i * "," * j * "," * NI])
)

plot_data_0_1$variable <- factor(
  plot_data_0_1$variable,
  levels = c(
    "cp_current_common_partner",
    "cp_continuous_cov",
    "cp_categorical_cov",
    "cp_number_interaction"
  )
)

plot_data_0_1 <- rbind(
  plot_data_0_1,
  plot_data_0_1[1, ]
)
plot_data_0_1$variable[length(plot_data_0_1$variable)] <- "cp_number_interaction"
plot_data_0_1$value[length(plot_data_0_1$value)] <- NA
N_const <- length(res_simulation[[1]]$popularity_1_0)
M_const <- mean(unlist(lapply(res_simulation, function(x) x$n_events_0_1)))
P_1 <- length(res_simulation[[1]]$est_0_1_list[[1]])
P_2 <- length(res_simulation[[1]]$est_1_0_list[[1]])

label_tmp <- TeX("$\\frac{{(P^{0 \\rightarrow 1})}^2}{M^{0 \\rightarrow 1}}$")

a <- ggplot(data = plot_data_0_1, aes(x = n_change, y = value, color = variable)) +
  theme_pubr(base_size = 17) +
  geom_vline(xintercept = 10, linetype = "dashed", linewidth = 1) +
  geom_smooth(na.rm = FALSE, show.legend = TRUE, se = FALSE) +
  geom_point() +
  ylim(c(0.89, 1)) +
  scale_color_manual(
    name = "Statistics",
    values = master_colors,
    labels = master_labels,
    breaks = names(master_colors),
    drop = FALSE
  ) +
  scale_x_continuous(
    name = "Number of steps of f (Q)",
    breaks = pretty(res_simulation[[1]]$grid_cells),
    sec.axis = sec_axis(
      transform = ~ (P_1 + N_const + .)^2 / (0.5 * M_const),
      breaks = round((pretty(res_simulation[[1]]$grid_cells) * 2 +
        length(res_simulation[[1]]$popularity_1_0) * 2)^2 /
        mean(unlist(lapply(res_simulation, function(x) x$n_events_0_1)))),
      name = label_tmp
    )
  ) +
  guides(color = guide_legend(override.aes = list(linetype = 1, shape = 16))) +
  xlab("Number of steps of f(Q)") +
  ylab("Coverage probability")
a



plot_data_1_0 <- melt.data.table(plot_data_1_0, id.vars = "n_change")
plot_data_1_0$variable <- factor(
  plot_data_1_0$variable,
  levels = c(
    "cp_current_common_partner",
    "cp_continuous_cov",
    "cp_categorical_cov",
    "cp_number_interaction"
  )
)
plot_data_1_0 <- rbind(
  plot_data_1_0,
  plot_data_1_0[1, ]
)
plot_data_1_0$variable[length(plot_data_1_0$variable)] <- "cp_current_common_partner"
plot_data_1_0$value[length(plot_data_1_0$value)] <- NA
label_tmp <- TeX("$\\frac{{(P^{1 \\rightarrow 0})}^2}{M^{1 \\rightarrow 0}}$")

M_const <- mean(unlist(lapply(res_simulation, function(x) x$n_events_1_0)))

b <- ggplot(data = plot_data_1_0, aes(x = n_change, y = value, color = variable)) +
  theme_pubr(base_size = 17) +
  geom_smooth(na.rm = FALSE, show.legend = TRUE, se = FALSE) +
  geom_point() +
  geom_vline(xintercept = 10, linetype = "dashed", linewidth = 1) +
  ylim(c(0.89, 1)) +
  scale_color_manual(
    name = "Statistics",
    values = master_colors,
    labels = master_labels,
    breaks = names(master_colors),
    drop = FALSE
  ) +
  scale_x_continuous(
    name = "Number of steps of f (Q)",
    breaks = pretty(res_simulation[[1]]$grid_cells),
    sec.axis = sec_axis(
      trans = ~ (P_2 + N_const + .)^2 / (0.5 * M_const),
      breaks = round((pretty(res_simulation[[1]]$grid_cells) * 2 +
        length(res_simulation[[1]]$popularity_1_0) * 2)^2 /
        mean(unlist(lapply(res_simulation, function(x) x$n_events_1_0)))),
      name = label_tmp
    )
  ) +
  guides(color = guide_legend(override.aes = list(linetype = 1, shape = 16))) +
  xlab("Number of steps of f(Q)") +
  ylab("Coverage probability")

final_plot <- (a | b) +
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = "(b) Coverage probability for varying steps of f(Q)",
    theme = theme(
      plot.caption = element_text(hjust = 0.5, size = 17)
    )
  ) &
  theme(legend.position = "bottom")

height_ratio <- c(1, 2)

ggsave(
  plot = plot_tmp / wrap_elements(full = final_plot) + plot_layout(
    heights = height_ratio
  ), filename = "Simulation/Plots/simulation_5.png",
  width = 10, height = 10
)
