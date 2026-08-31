
rm(list=ls())
source("Simulation/helper.R")
library(stringr)
library(readr)
library(stringi)
library(MASS)
library(redeem)
library(Hmisc)
library(patchwork)
library(parallel)
library(grid)
library(gridExtra)
library(ggpubr)
load(file = "Simulation/Results/setting_simulation_2.RData")

res_simulation <- readRDS(file = "Simulation/Results/simulation_2.RDS")
coef_0_1 = res_simulation[[1]]$coef_0_1
coef_1_0 = res_simulation[[1]]$coef_1_0


error_res <- get_error(res_simulation,  coef_0_1 = res_simulation[[1]]$coef_0_1,
                       popularity_0_1= res_simulation[[1]]$popularity_0_1,
                       baseline_0_1_gt = res_simulation[[1]]$baseline_0_1_gt,
                       coef_1_0 = res_simulation[[1]]$coef_1_0,
                       popularity_1_0 = res_simulation[[1]]$popularity_1_0,
                       baseline_1_0_gt = res_simulation[[1]]$baseline_1_0_gt)


se_0_1 <- do.call("rbind",lapply(res_simulation, function(x) sqrt(diag(x$covariance_0_1))))
se_1_0 <- do.call("rbind",lapply(res_simulation, function(x) sqrt(diag(x$covariance_1_0))))
est_0_1 <- t(do.call("cbind",lapply(res_simulation, function(x) x$coef_0_1_core)))
est_1_0 <- t(do.call("cbind",lapply(res_simulation, function(x) x$coef_1_0_core)))

plot_list <- list()

names_tmp <-  c("Current Common Part.", "Continuous Cov.", "Categorical Cov.",
                "Number Interaction", "Continuous Cov.", "Categorical Cov.")

size_tmp <- 20
for(i in seq_len(ncol(est_0_1))){
  x <- (est_0_1[,i]- res_simulation[[1]]$coef_0_1[i])/se_0_1[,i]
  data_tmp <- qqnorm(x)
  probs <- c(0.25,0.75)
  y <- as.vector(quantile(x, probs, names = FALSE, type = 7,
                          na.rm = TRUE))
  x <- qnorm(probs)
  slope <- diff(y)/diff(x)
  int <- y[[1L]] - slope * x[[1L]]

  data_tmp <- as.data.frame(data_tmp)

  plot_list[[i]] <- ggplot(data = data_tmp, aes(x = x, y =y)) +
    geom_point(alpha = 0.1) +
    theme_pubr(base_size = 17) +
    theme(axis.title.y = element_text(size = size_tmp)) +
    ylab("Sample quantities") +
    xlab(" ") +
    geom_abline(intercept = int,slope = slope) +
    ggtitle(names_tmp[i]) +
    theme(plot.title = element_text(hjust = 0.5))


}

for(i in seq_len(ncol(est_1_0))){
  x <- (est_1_0[,i]- res_simulation[[1]]$coef_1_0[i])/se_1_0[,i]
  data_tmp <- qqnorm(x)
  probs <- c(0.25,0.75)
  y <- as.vector(quantile(x, probs, names = FALSE, type = 7,
                          na.rm = TRUE))
  x <- qnorm(probs)
  slope <- diff(y)/diff(x)
  int <- y[[1L]] - slope * x[[1L]]

  data_tmp <- as.data.frame(data_tmp)
  plot_list[[i+ ncol(est_0_1)]] <- ggplot(data = data_tmp, aes(x = x, y =y)) +
    geom_point(alpha = 0.1) +
    theme_pubr(base_size = 17) +
    theme(axis.title.y = element_text(size = size_tmp)) +
    ylab("Sample quantities") +
    xlab(" ") +
    geom_abline(intercept = int,slope = slope) +
    ggtitle(names_tmp[i + ncol(est_0_1)]) +
    theme(plot.title = element_text(hjust = 0.5))

}

grid1 <- grid.arrange(grobs = list(plot_list[[1]], plot_list[[2]], plot_list[[3]]),
                      bottom = textGrob("Theoretical Quantities", gp=gpar(fontsize=size_tmp), hjust = 0.3,vjust = -1),
             top = textGrob("Incidence Model", gp=gpar(fontsize=25)), ncol=3)
grid2 <- grid.arrange(grobs = list(plot_list[[4]], plot_list[[5]], plot_list[[6]]),
                      bottom = textGrob("Theoretical Quantities", gp=gpar(fontsize=size_tmp), hjust = 0.3,vjust = -1),
             top = textGrob("Duration Model", gp=gpar(fontsize=25)), ncol=3)

Cairo::CairoPNG("Simulation/Plots/normality_estimators.png",
                width = 720*4,height = 500*4, res = 120*2)
grid.newpage()
grid.draw(rbind(grid1, grid2, size = "last"))
dev.off()


table_res <- round(rbind(cbind(error_res$ave_0_1 - coef_0_1, error_res$error_rmse_core_0_1,error_res$cp_0_1),
                   cbind(error_res$ave_1_0 - coef_1_0, error_res$error_rmse_core_1_0,error_res$cp_1_0)),3)

colnames(table_res) <- c("Bias", "RMSE", "CP")
rownames(table_res) <- c("s_{i,j,CCP}", "s_{i,j,\bm{x}}", "s_{i,j,\bm{w}}",
                         "s_{i,j,NI}", "s_{i,j,\bm{x}}", "s_{i,j,\bm{w}}")
# Info for table 1 in the simulation study
latex_res <- cbind(cbind(rownames(table_res)[1:3], table_res[1:3,]),
      cbind(rownames(table_res)[4:6], table_res[4:6,]))
latex_res <- latex_res[,-1]

new_res <- latex_res
i = 1
for(i in 1:ncol(latex_res)){
  if(i != 4){
    new_res[,i] <- sprintf("%.3f", as.numeric(latex_res[,i]))
    new_res[,i] <- gsub(pattern = "0\\.",replacement = ".",new_res[,i])
  }

}

new_res <- gsub(pattern = "-",replacement = "$-$",new_res)
new_res[,4] <- paste0("$",new_res[,4],"$")
rownames(new_res) <- paste0("$",rownames(new_res),"$")
rownames(new_res) <- gsub(pattern = "\bm",replacement = "\\\\bm",rownames(new_res))
new_res[,4] <- gsub(pattern = "\bm",replacement = "\\\\bm",new_res[,4])

latex(object = new_res,n.cgroup = c(3,4),
      cgroup = c("Incidence ($\\alpha^{0\\rightarrow 1}$)",
                 "Duration ($\\alpha^{1\\rightarrow 0}$)"),
      file = "Simulation/Results/table_simulation_1.tex",
      label = "tbl:simulation",
      title = "",
      first.hline.double= F,
      rowlabel.just = c("l"),
      caption = paste("For each effect, we report the Bias ($\\hat{\\bm{\\alpha}}-\\bm{\\alpha}^\\star$), RMSE (Root-Mean-Squared Error), and CP (Coverage Probability). The dependence of the sufficient statistics on the history $\\mathscr{H}_t$ is suppressed for notational clarity.."),
      cgroup.just="c")

degree_0_1 = res_simulation[[1]]$popularity_0_1
plot_data <- data.frame(degree_0_1 = degree_0_1,
                        error_rmse_degree_0_1 =error_res$error_inf_degree_0_1)

Cairo::CairoPNG("Simulation/Plots/error_pop_0_1.png",
                width = 692*2,height = 692*2, res = 120*2)
ggplot(plot_data, aes(x = degree_0_1, y = error_rmse_degree_0_1)) +
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)),
                  bins = 10,n = 500) +
  scale_fill_distiller("",palette = "Blues", direction = 1) +
  labs(
    x = expression(beta^{0 %->% 1}),
    y = expression(RMSE * " " * bgroup("(", hat(beta)^{0 %->% 1}, ")"))
  ) +
  theme_pubr()

dev.off()

degree_1_0 = res_simulation[[1]]$popularity_1_0
plot_data <- data.frame(degree_0_1 = degree_1_0,
                        error_rmse_degree_0_1 =error_res$error_rmse_degree_1_0)
library(ggpubr)
Cairo::CairoPNG("Simulation/Plots/error_pop_1_0.png",
                width = 692*2,height = 692*2, res = 120*2)
ggplot(plot_data, aes(x = degree_0_1, y = error_rmse_degree_0_1)) +
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)),
                  bins = 10,n = 500) +
  scale_fill_distiller("",palette = "Blues", direction = 1) +
  labs(
    x = expression(beta^{1 %->% 0}),
    y = expression(RMSE * " " * bgroup("(", hat(beta)^{1 %->% 0}, ")"))
  ) +
  theme_pubr()
dev.off()

plot_data <- data.frame(error_rmse_degree_0_1 = error_res$error_rmse_degree_0_1,
                        error_rmse_degree_1_0 =error_res$error_rmse_degree_1_0)
library(ggpubr)
Cairo::CairoPNG("Simulation/Plots/error_pop_both.png",
                width = 692*2,height = 692*2, res = 120*2)
ggplot(plot_data, aes(x = error_rmse_degree_0_1, y = error_rmse_degree_1_0)) +
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)),
                  bins = 10,n = 500) +
  scale_fill_distiller("",palette = "Blues", direction = 1) +
  labs(
    x = expression(beta^{1 %->% 0}),
    y = expression(RMSE * " " * bgroup("(", hat(beta)^{1 %->% 0}, ")"))
  ) +
  theme_pubr()
dev.off()

