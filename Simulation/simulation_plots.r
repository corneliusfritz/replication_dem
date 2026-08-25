
rm(list=ls())
source("Simulation/helper.R")
library(stringr)
library(readr)
library(stringi)
library(MASS)
library(redeem)
library(parallel)
# library(lineprof)
library(peakRAM)
library(ggpubr)
library(ggpubr)
library(dplyr)
library(cowplot)
library(patchwork)
library(scales)
set.seed(123456)

label_clean_numbers <- function(x) {
  sapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    if (abs(v) >= 1000) return(comma(v, accuracy = 1))
    if (v %% 1 == 0) return(formatC(v, format = "f", digits = 0))
    s <- number(v, accuracy = 0.1, big.mark = ",", trim = TRUE)
    sub("^(-?)0\\.", "\\1.", s)
  })
}


n_nodes_poss <- c(50, 75, 100, 125, 150, 175)

simultaneous_interactions <- TRUE
it_max <- 1500
accelerated <- FALSE
estimate_newton <- FALSE
K <- 50
verbose <- FALSE


res <- list()
res_peak_ram <- list()
corr_data_list <- list()

errors <- list()
i <- 1

for(n_nodes in n_nodes_poss){
  cat("Loading ",n_nodes,"\n")
  tmp <- readRDS(paste0("Simulation/Results/simulation_3_",n_nodes, ".RDS"))
  # corr_data_list[[i]] <- unlist(lapply(tmp[[1]],FUN = function(x) {
  #   unidentifiable_0_1 <- x$unidentifiable_0_1
  #   if(length(unidentifiable_0_1)>0){
  #     x$est_mm_0_1[unidentifiable_0_1+3]  <- NA
  #     x$est_nr_0_1[unidentifiable_0_1+3]  <- NA
  #   }
  #   unidentifiable_1_0 <- x$unidentifiable_1_0
  #   if(length(unidentifiable_1_0)>0){
  #     x$est_mm_1_0[unidentifiable_1_0+3]  <- NA
  #     x$est_nr_1_0[unidentifiable_1_0+3]  <- NA
  #   }
  #   return(cor(c(x$est_mm_0_1,x$est_mm_1_0),c(x$est_nr_0_1,x$est_nr_1_0),use = "complete.obs"))
  # }))

  times_mm <- unlist(lapply(tmp[[1]],FUN = function(x) x$time_mm))
  times_nr <- unlist(lapply(tmp[[1]],FUN = function(x) x$time_nr))
  allocation_mm <- unlist(lapply(tmp[[1]],FUN = function(x) x$peak_ram_mm))
  allocation_nr <- unlist(lapply(tmp[[1]],FUN = function(x) x$peak_ram_nr))
  res[[i]] <- melt.data.table(data.table("Three-Step Estimation" = times_mm,
                                         "Newton Raphson" = times_nr,
                                         n_nodes = n_nodes),id.vars = "n_nodes")
  res_peak_ram[[i]]<- melt.data.table(data.table("Three-Step Estimation" = allocation_mm,
                                                 "Newton Raphson" = allocation_nr,
                                                 n_nodes = n_nodes),id.vars = "n_nodes")
  i <- i+1
}

res_temporal <- do.call(rbind, res)
res_peak_ram <- do.call(rbind, res_peak_ram)
# plots for simulation study 2
n_nodes_poss <- c(50, 100, 300, 500, 700, 1000)

simultaneous_interactions <- TRUE
it_max <- 1500
accelerated <- FALSE
estimate_newton <- FALSE
K <- 50
verbose <- FALSE

res <- list()
errors <- list()
n_events <- list()
i <- 1
n_nodes <- n_nodes_poss[5]
n_nodes <- n_nodes_poss[1]
for(n_nodes in n_nodes_poss){
  cat("Loading ",n_nodes,"\n")
  res[[i]] <- readRDS(paste0("Simulation/Results/simulation_1_",n_nodes, ".RDS"))
  # res1 <- readRDS(paste0("Simulation/Results/simulation_1_50.RDS"))
  # res2 <- readRDS(paste0("Simulation/Results/simulation_1_1000.RDS"))
  # object.size(res1[[1]][[1]])
  # object.size(res[[1]][[1]])
  # debugonce(get_error)
  errors[[i]] <- get_error(res[[i]][[1]],coef_0_1 = res[[i]]$coef_0_1,
                           popularity_0_1 = res[[i]]$popularity_0_1,
                           baseline_0_1_gt = res[[i]]$baseline_0_1_gt,
                           coef_1_0 = res[[i]]$coef_1_0,
                           popularity_1_0 = res[[i]]$popularity_1_0,
                           baseline_1_0_gt = res[[i]]$baseline_1_0_gt)

  n_events[[i]] <- unlist(lapply(res[[i]][[1]], function(x) x$number_events))

  cat(mean(errors[[i]]$error_inf_0_1),"\n")
  cat(mean(errors[[i]]$error_inf_1_0),"\n")
  cat(mean(errors[[i]]$error_inf_core_0_1),"\n")
  cat(mean(errors[[i]]$error_inf_core_1_0),"\n")
  cat(mean(errors[[i]]$error_inf_time_0_1),"\n")
  cat(mean(errors[[i]]$error_inf_time_1_0),"\n")


  # res[[1]]$
  i <- i+1
}


plot_data <- data.table(error = c(errors[[1]]$error_rmse_alt_degree_0_1,
                                  errors[[2]]$error_rmse_alt_degree_0_1,
                                  errors[[3]]$error_rmse_alt_degree_0_1,
                                  errors[[4]]$error_rmse_alt_degree_0_1,
                                  errors[[5]]$error_rmse_alt_degree_0_1,
                                  errors[[6]]$error_rmse_alt_degree_0_1,
                                  errors[[1]]$error_rmse_alt_degree_1_0,
                                  errors[[2]]$error_rmse_alt_degree_1_0,
                                  errors[[3]]$error_rmse_alt_degree_1_0,
                                  errors[[4]]$error_rmse_alt_degree_1_0,
                                  errors[[5]]$error_rmse_alt_degree_1_0,
                                  errors[[6]]$error_rmse_alt_degree_1_0,
                                  errors[[1]]$error_rmse_alt_core_0_1,
                                  errors[[2]]$error_rmse_alt_core_0_1,
                                  errors[[3]]$error_rmse_alt_core_0_1,
                                  errors[[4]]$error_rmse_alt_core_0_1,
                                  errors[[5]]$error_rmse_alt_core_0_1,
                                  errors[[6]]$error_rmse_alt_core_0_1,
                                  errors[[1]]$error_rmse_alt_core_1_0,
                                  errors[[2]]$error_rmse_alt_core_1_0,
                                  errors[[3]]$error_rmse_alt_core_1_0,
                                  errors[[4]]$error_rmse_alt_core_1_0,
                                  errors[[5]]$error_rmse_alt_core_1_0,
                                  errors[[6]]$error_rmse_alt_core_1_0,
                                  errors[[1]]$error_rmse_alt_time_0_1,
                                  errors[[2]]$error_rmse_alt_time_0_1,
                                  errors[[3]]$error_rmse_alt_time_0_1,
                                  errors[[4]]$error_rmse_alt_time_0_1,
                                  errors[[5]]$error_rmse_alt_time_0_1,
                                  errors[[6]]$error_rmse_alt_time_0_1,
                                  errors[[1]]$error_rmse_alt_time_1_0,
                                  errors[[2]]$error_rmse_alt_time_1_0,
                                  errors[[3]]$error_rmse_alt_time_1_0,
                                  errors[[4]]$error_rmse_alt_time_1_0,
                                  errors[[5]]$error_rmse_alt_time_1_0,
                                  errors[[6]]$error_rmse_alt_time_1_0),
                        n_nodes = rep(rep(n_nodes_poss, each = length(errors[[5]]$error_inf_1_0)), times = 6),
                        model = rep(rep(c("Incidence", "Duration"), each = 6* length(errors[[5]]$error_inf_1_0)), times = 3),
                        type = rep(c("Popularity", "Non-Popularity","Temporal"), each = length(errors[[5]]$error_inf_1_0)*12))


size <- 22
means <- res_temporal[,.(mean_value = mean(value)), by = .(n_nodes,variable)]
p_1 <- ggplot(res_temporal, aes(y = value, color = variable, x = factor(n_nodes)))+
  geom_boxplot(position = position_dodge(width = 0),show.legend = T) +
  geom_line(data = means, aes(x = factor(n_nodes), y = mean_value, group = variable),
            show.legend = FALSE,lty =1) +
  labs(caption = "(b) Simulation Study 3: Time Comparison") +
  scale_color_manual("Algorithm  ",values = c("black", "grey"))+
  xlab("Number of Actors (N)")+
  ylab("Time in Seconds")+
  theme_pubr(base_size = size)+
  scale_y_continuous(labels = label_clean_numbers) +
  theme(plot.margin = margin(0, 0, 0, 0),
        plot.caption = element_text(size = size))

p_1

means <- res_peak_ram[,.(mean_value = mean(value)), by = .(n_nodes,variable)]

p_2 <- ggplot(res_peak_ram, aes(y = value, color = variable, x = factor(n_nodes)))+
  geom_boxplot(position = position_dodge(width = 0),show.legend = T) +
  geom_line(data = means, aes(x = factor(n_nodes), y = mean_value, group = variable),
            show.legend = FALSE,lty =1) +
  scale_color_manual("Algorithm  ",values = c("black", "grey"))+
  labs(caption = "(c) Simulation Study 3: Memory Comparison") +
  xlab("Number of Actors (N)")+
  ylab("Peak Memory in MB")+
  theme_pubr(base_size = size) +
  scale_y_continuous(labels = label_clean_numbers) +
  theme(plot.margin = margin(0, 0, 0, 0),
        plot.caption = element_text(size = size))

p_2
# plot_data$by <- paste(plot_data$model, plot_data$t)


p_0 <- ggplot(plot_data, aes(y = error, x = factor(n_nodes), color = type))+
  geom_boxplot(position = position_dodge2(width = 0.75, preserve = "single"), alpha = 0.06) +
  xlab("Number of Actors (N)")+
  ylab("RMSE")+
  theme_pubr(base_size = size) +
  labs(caption = "(a) Simulation Study 2: Estimation Error under Increasing Actors") +
  scale_y_continuous(labels = label_clean_numbers) +
  scale_x_discrete(labels = label_clean_numbers(unique(plot_data$n_nodes))) +
  scale_color_manual("Estimator  ",values = c("#0072B2", "#E69F00", "#D55E00"),
                     labels = c(expression(hat(alpha)),
                                expression(hat(beta)),
                                expression(hat(gamma)))) +
  theme(axis.title.y = element_text(size = size),
        plot.caption = element_text(size = size, hjust = 0.5),
        # legend.text = element_text(size = size),
        legend.position      = c(0.6, 1.05),
        legend.justification = c(1, 1),
        legend.direction = "horizontal")

p_3 <- ggarrange(
  p_1+ theme(legend.text = element_text(margin = margin(r = 10, l = 10))) +
  guides(color = guide_legend(
    keyheight = unit(1.2, "cm")
  )),
  p_2+ theme(legend.text = element_text(margin = margin(r = 10, l = 10)))+
    guides(color = guide_legend(
      keyheight = unit(1.2, "cm")
    )), ncol = 2,
  common.legend = TRUE, legend = "top",
  align = "v", heights = c(1, 1)
)

ggsave(plot = plot_grid(
  p_0 +theme(plot.margin = margin(t = 10, r = 5.5, b = 5.5, l = 25, unit = "pt")), p_3,
  ncol = 1,rel_heights = c(1, 1)
), "Simulation/Plots/simulation_final.png",
       width = 14,height = 13)


