library(data.table)
library(redeem)
library(stringr)
library(ggpubr)
library(ggExtra)
library(scales)
library(latex2exp)
library(Hmisc)
library(patchwork)

calls = fread("Application/Data/calls.csv")
actor_data = fread("Application/Data/actor_data.csv")

label_clean_numbers <- function(x) {
  sapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    if (abs(v) >= 1000) return(comma(v, accuracy = 1))
    if (v %% 1 == 0) return(formatC(v, format = "f", digits = 0))
    s <- number(v, accuracy = 0.1, big.mark = ",", trim = TRUE)
    sub("^(-?)0\\.", "\\1.", s)
  })
}
# sort(unique(calls$caller))
# sort(unique(actor_data$unique_actors))
calls$caller <- match(calls$caller, actor_data$unique_actors)
calls$callee <- match(calls$callee, actor_data$unique_actors)
# Only look at calls between users of the study and with a duration larger than 0
calls <- calls[!is.na(caller) &!is.na(callee) & caller != callee  & duration >0]
names(calls) <- c("timestamp", "from", "to", "duration")
calls$time <- calls$timestamp - min(calls$timestamp)
# Make sure that i<j holds for the undirected durational events
calls$tmp <- calls$to<calls$from
tmp <- calls$to
calls$to[calls$tmp] <- calls$from[calls$tmp]
calls$from[calls$tmp] <- tmp[calls$tmp]
calls$tmp <- NULL
calls$from_to <- paste(calls$from, calls$to, sep = "_")
calls <- calls[,.(time,duration, from, to,from_to)]

tmp_fun <- function(x) {
  start <- x$time
  end <- x$time + x$duration
  res <- data.frame(rbind(cbind(x$from, x$to, start, 1),
        cbind(x$from, x$to, end, 0)))
  colnames(res) <- c("from", "to", "time", "type")
  res <- res[order(res$time),]
  if(length(which(diff(res$type) == 0))>0){
    number_errors <- length(which(diff(res$type) == 0))/2

    errors <- which(diff(res$type) == 0)
    exclude <- c(errors, errors+1)

    exclude <- c()
    for(i in seq(from = 1, to = length(errors), by = 2)){
      exclude <- c(exclude, errors[i]:(errors[i+1]+1))
    }

    res <- rbind(res[errors[seq(1,length(errors),2)],],
          res[errors[seq(2,length(errors),2)] +1,],
          res[-exclude,])

    res <- res[order(res$time),]
    #
    #
    #
    # return(NA)
  }
  return(res)
}

res <- calls[,.(mat = list(tmp_fun(.SD))), by = from_to]
# res$from_to[which(is.na(res$mat))]

res <- data.table(do.call(what = rbind,res$mat))

events <- res[,.(time, from, to, type)]

events$time <- as.numeric(events$time)
events$time <- events$time - min(events$time)
# events$time = 1:length(events$time)
colnames(actor_data) <- c("id","id_old")
proximity <- fread("Application/Data/events.csv")

gender <- fread("Application/Data/genders.csv")
actor_data$female <- gender$female[match(actor_data$id_old,gender$`# user`)]

sms <-  fread("Application/Data/sms.csv")

fb_friends <-  fread("Application/Data/fb_friends.csv")
names(fb_friends) <- c("from", "to")
fb_friends <- fb_friends[!from == to]
fb_friends <- fb_friends[!duplicated(fb_friends)]
fb_friends$included <- (fb_friends$from %in% actor_data$id_old) & (fb_friends$to %in% actor_data$id_old)
fb_friends <- fb_friends[included == T]
fb_friends$from_id <- match(fb_friends$from,actor_data$id_old)
fb_friends$to_id <- match(fb_friends$to,actor_data$id_old)
friendship_mat <- matrix(0,nrow = nrow(actor_data),ncol = nrow(actor_data))
friendship_mat[cbind(fb_friends$from_id, fb_friends$to_id)] = 1
friendship_mat[cbind(fb_friends$to_id, fb_friends$from_id)] = 1

both_female <- outer(actor_data$female, actor_data$female, FUN = function(x,y){x==y})
one_female <- outer(actor_data$female, 1- actor_data$female, FUN = function(x,y){x==y})
actor_data$exclude <- !actor_data$id %in% unique(as.numeric(cbind(events$from,events$to)))
# !actor_data$id %in% unique(as.numeric(cbind(proximity$from,proximity$to)))
actor_data$exclude[is.na(actor_data$female)] <- TRUE
actor_data$exclude[which(rowSums(friendship_mat) ==0)] <- TRUE

for(i in 1:10){
  events_cut <- events[events$from %in% actor_data$id[!actor_data$exclude] &
                        events$to %in% actor_data$id[!actor_data$exclude]]

  actor_data$exclude[which(!actor_data$id %in% unique(as.numeric(cbind(events_cut$from,events_cut$to))))] <- TRUE

  proximity_cut <- proximity[proximity$from %in% actor_data$id[!actor_data$exclude] &
                              proximity$to %in% actor_data$id[!actor_data$exclude]]

  actor_data$exclude[which(!actor_data$id %in% unique(as.numeric(cbind(proximity_cut$from,proximity_cut$to))))] <- TRUE
  cat(sum(!actor_data$exclude),"\n")
}



friendship_mat_cut <- friendship_mat[!actor_data$exclude, !actor_data$exclude]
both_female_cut <- both_female[!actor_data$exclude, !actor_data$exclude]
one_female_cut <- one_female[!actor_data$exclude, !actor_data$exclude]

actor_data$id_included <- NA
actor_data$id_included[!actor_data$exclude] <- 1:sum(!actor_data$exclude)
proximity_cut$from <- actor_data$id_included[match(proximity_cut$from,actor_data$id)]
proximity_cut$to <- actor_data$id_included[match(proximity_cut$to,actor_data$id)]

events_cut$from <- actor_data$id_included[match(events_cut$from,actor_data$id)]
events_cut$to <- actor_data$id_included[match(events_cut$to,actor_data$id)]

sort(unique(c(events_cut$from, events_cut$to)))

# daily fixed effect
time_changepoints <- seq(from = 60*60*24, by = 60*60*24, to = max(proximity$time))
time_changepoints <- time_changepoints[-length(time_changepoints)]
# time_changepoints <- unique(proximity$time)
# time_changepoints <- time_changepoints[-length(time_changepoints)]
labels_changepoints <- c(1:(length(time_changepoints)))
events_cut <- events_cut[order(time)]
events_cut <- as.matrix(events_cut)
n_nodes <- nrow(friendship_mat_cut)

options(expressions = 500000)

models <- readRDS("Application/Results/denmark_models_call.rds")
residuals <- get_residuals(models)

plot(residuals$surv, residuals$theoretical)
models_proximity <- readRDS("Application/Results/denmark_models_proximity.rds")
residuals_proximity <- get_residuals(models_proximity)
models_alt <- readRDS("Application/Results/denmark_models_call_alt.rds")
summary(models_alt)


models_proximity_alt <- readRDS("Application/Results/denmark_models_proximity_alt.rds")
#  Degree structure -------


event_start <- events_cut[events_cut[,4]==1,]
proximity_start <- as.matrix(proximity_cut)[as.matrix(proximity_cut)[,4]==1,]

dens_call <- density(table(as.numeric(event_start[,c(2,3)])),
                     from = 0, n = 1000)
dens_call <- data.frame(x = dens_call$x,
                        y = dens_call$y)

dens_proximity <- density(table(as.numeric(proximity_start[,c(2,3)])),
                          from = 0, n = 1000)
dens_proximity <- data.frame(x = dens_proximity$x,
                             y = dens_proximity$y)

my_theme <- theme_pubr(20) +
  theme(plot.title = element_text(hjust = 0.5), legend.key.size = unit(1.5, "cm"))

a <- ggplot(data = dens_call, aes(x = x, y = y)) +
  geom_line() + my_theme + xlab("Number of cumulative interactions") +
  ylab("Density")

b <- ggplot(data = dens_proximity, aes(x = x, y = y)) +
  geom_line() + my_theme + xlab("Number of cumulative interactions") +
  ylab("Density")

Cairo::CairoPNG("Application/Plots/degrees_a.png",
                width = 692*2,height = 692*2, res = 120*2)
a
dev.off()

Cairo::CairoPNG("Application/Plots/degrees_b.png",
                width = 692*2,height = 692*2, res = 120*2)
b
dev.off()
a <- ggplot(data = dens_call, aes(x = x, y = y)) +
  geom_line() + my_theme + xlab("Number of cumulative interactions") +
  ylab("Density") + ggtitle("Digital (Call)")

b <- ggplot(data = dens_proximity, aes(x = x, y = y)) +
  geom_line() + my_theme + xlab("Number of cumulative interactions") +
  ylab("Density") + ggtitle("Physical (Co-location)")

Cairo::CairoPNG("Application/Plots/degrees.png",
                width = 692*4,height = 692*2, res = 120*2)

b/a
dev.off()


#  Baseline call -------

time <- as.numeric(gsub(pattern = "time_", "", names(models$model_0_1$est_time)))
time <- time/60/60/24 #transform time from seconds to days

step_data_1 <- data.frame(
  x = c(0,time),
  y = c(0,models$model_0_1$est_time)
)

smoothed <- loess(step_data_1$y~ step_data_1$x)
smooth_data_2 <- data.frame(
  x = step_data_1$x,
  y = smoothed$fitted
)

step_data_2 <- data.frame(
  x = step_data_1$x,
  y = smoothed$residuals
)


a <- ggplot(step_data_1, aes(x = x, y = y)) +
  geom_step(direction = "hv", colour = "gray")  +
  geom_smooth(se = FALSE,colour = "black")+
  theme_pubr(base_size = 20) +
  ylim(c(range(step_data_1$y)))+
  ylab(expression(f(t,hat(gamma)^{0 %->% 1}))) +
  xlab("t in Days")
a

c <- ggplot(step_data_2, aes(x = x, y = y)) +
  geom_step(colour = "gray")  +
  # ylim(c(range(step_data_1$y)))+
  theme_pubr(base_size = 20) +
  ylab("Residual of trend component") +
  xlab("t in days")

Cairo::CairoPNG("Application/Plots/baseline_call_0_1.png",
                width = 1011*2,height = 692*4, res = 120*2)


a/c
dev.off()

time <- as.numeric(gsub(pattern = "time_", "", names(models$model_1_0$est_time)))
time <- time/60/60/24 #transform time from seconds to days

step_data_1 <- data.frame(
  x = c(0,time),
  y = c(0,models$model_1_0$est_time)
)

smoothed <- loess(step_data_1$y~ step_data_1$x)
smooth_data_2 <- data.frame(
  x = step_data_1$x,
  y = smoothed$fitted
)

step_data_2 <- data.frame(
  x = step_data_1$x,
  y = smoothed$residuals
)


a <- ggplot(step_data_1, aes(x = x, y = y)) +
  geom_step(direction = "hv", colour = "gray")  +
  geom_smooth(se = FALSE,colour = "black")+
  theme_pubr(base_size = 20) +
  ylim(c(range(step_data_1$y)))+
  ylab(expression(f(t,hat(gamma)^{1 %->% 0}))) +
  xlab("t in days")

c <- ggplot(step_data_2, aes(x = x, y = y)) +
  geom_step(colour = "gray")  +
  # ylim(c(range(step_data_1$y)))+
  theme_pubr(base_size = 20) +
  ylab("Residual of trend component") +
  xlab("time in days (t)")



Cairo::CairoPNG("Application/Plots/baseline_call_1_0.png",
                width = 1011*2,height = 692*4, res = 120*2)

a/c
dev.off()

#  Baseline Physical (Co-location) -------
time <- as.numeric(gsub(pattern = "time_", "", names(models_proximity$model_0_1$est_time)))
time <- time/60/60/24 #transform time from seconds to days

step_data_1 <- data.frame(
  x = c(0,time),
  y = c(0,models_proximity$model_0_1$est_time)
)

smoothed <- loess(step_data_1$y~ step_data_1$x)
smooth_data_2 <- data.frame(
  x = step_data_1$x,
  y = smoothed$fitted
)

step_data_2 <- data.frame(
  x = step_data_1$x,
  y = smoothed$residuals
)


a <- ggplot(step_data_1, aes(x = x, y = y)) +
  geom_step(direction = "hv", colour = "gray")  +
  geom_smooth(se = FALSE,colour = "black")+
  theme_pubr(base_size = 20) +
  ylim(c(range(step_data_1$y)))+
  ylab(expression(f(t,hat(gamma)^{0 %->% 1}))) +
  xlab("time in days (t)")

c <- ggplot(step_data_2, aes(x = x, y = y)) +
  geom_step(colour = "gray")  +
  # ylim(c(range(step_data_1$y)))+
  theme_pubr(base_size = 20) +
  ylab("Residual of trend component") +
  xlab("time in days (t)")

Cairo::CairoPNG("Application/Plots/baseline_proximity_0_1.png",
                width = 1011*2,height = 692*4, res = 120*2)

a/c
dev.off()

time <- as.numeric(gsub(pattern = "time_", "", names(models_proximity$model_1_0$est_time)))
time <- time/60/60/24 #transform time from seconds to days

step_data_1 <- data.frame(
  x = c(0,time),
  y = c(0,models_proximity$model_1_0$est_time)
)

smoothed <- loess(step_data_1$y~ step_data_1$x)
smooth_data_2 <- data.frame(
  x = step_data_1$x,
  y = smoothed$fitted
)

step_data_2 <- data.frame(
  x = step_data_1$x,
  y = smoothed$residuals
)


a <- ggplot(step_data_1, aes(x = x, y = y)) +
  geom_step(direction = "hv", colour = "gray")  +
  geom_smooth(se = FALSE,colour = "black")+
  theme_pubr(base_size = 20) +
  ylim(c(range(step_data_1$y)))+
  ylab(expression(f(t,hat(gamma)^{1 %->% 0}))) +
  xlab("time in days (t)")


c <- ggplot(step_data_2, aes(x = x, y = y)) +
  geom_step(colour = "gray")  +
  theme_pubr(base_size = 20) +
  ylab("Residual of trend component") +
  xlab("time in days (t)")



Cairo::CairoPNG("Application/Plots/baseline_proximity_1_0.png",
                width = 1011*2,height = 692*4, res = 120*2)
a/c
dev.off()



exo_changes = seq(from = 60*60*24, by = 60*60*24, to = max(proximity$time))
exo_labels = c(1,1:(length(exo_changes)))

n_nodes = nrow(friendship_mat_cut)


# Popularity Estimates ----
n_nodes = nrow(friendship_mat_cut)
training_start <- 0
proximity_train <- proximity_cut[proximity_cut$time>training_start,]
identifiable <- unique(sort(c(proximity_train$to,
                              proximity_train$from)))
effects_pop <- data.frame(coef_0_1 = models_proximity$model_0_1$est_degree[identifiable],
                          coef_1_0 = models_proximity$model_1_0$est_degree[identifiable])

p <- ggplot(effects_pop, aes(x = coef_0_1, y = coef_1_0)) +
  geom_density2d(color = "black",alpha = 0.7) +
  stat_density_2d(geom = "polygon", alpha = 0.1) +
  geom_point(alpha = 0) +
  theme_pubr(base_size = 30) +
  xlab(expression(hat(beta)^{0 %->% 1})) +
  ylab(expression(hat(beta)^{1 %->% 0})) +
  xlim(c(-8.5,-7.2)) +
  ylim(c(-4,-3))

Cairo::CairoPNG("Application/Plots/effeect_pop_proximity.png",
                width = 1011*2,height = 692*2, res = 120*2)
ggMarginal(
  p,
  type = "density",    # Add density plots on the margins
  fill = "black", alpha =0.2
)

dev.off()


effects_pop <- data.frame(coef_0_1 = models$model_0_1$est_degree[identifiable],
                          coef_1_0 = models$model_1_0$est_degree[identifiable])

p <- ggplot(effects_pop, aes(x = coef_0_1, y = coef_1_0)) +
  geom_density2d(color = "black",alpha = 0.7) +
  stat_density_2d(geom = "polygon", alpha = 0.1) +
  geom_point(alpha = 0) +
  theme_pubr(base_size = 30) +
  xlab(expression(hat(beta)^{0 %->% 1})) +
  ylab(expression(hat(beta)^{1 %->% 0})) +
  xlim(c(-12.5,-9)) +
  ylim(c(-2.5,4))


Cairo::CairoPNG("Application/Plots/effeect_pop_call.png",
                width = 1011*2,height = 692*2, res = 120*2)
ggMarginal(
  p,
  type = "density",    # Add density plots on the margins
  fill = "black", alpha =0.2
)
dev.off()

effects_pop <- data.frame(call = models$model_0_1$est_degree[identifiable],
                          proximity = models_proximity$model_0_1$est_degree[identifiable])

p <- ggplot(effects_pop, aes(x = call, y = proximity)) +
  geom_density2d(color = "black",alpha = 0.7) +
  stat_density_2d(geom = "polygon", alpha = 0.1) +
  geom_point(alpha = 0.3) +
  theme_pubr(base_size = 30) +
  xlab("Digital (Call)") +
  ylab("Physical (Co-location)")

p

Cairo::CairoPNG("Application/Plots/effect_pop_call_proxy_0_1.png",
                width = 1011*2,height = 692*2, res = 120*2)
p
dev.off()

effects_pop <- data.frame(call = models$model_1_0$est_degree[identifiable],
                          proximity = models_proximity$model_1_0$est_degree[identifiable])

p <- ggplot(effects_pop, aes(x = call, y = proximity)) +
  geom_density2d(color = "black",alpha = 0.7) +
  stat_density_2d(geom = "polygon", alpha = 0.1) +
  geom_point(alpha = 0.3) +
  theme_pubr(base_size = 30) +
  xlab("Digital (Call)") +
  ylab("Physical (Co-location)")

Cairo::CairoPNG("Application/Plots/effect_pop_call_proxy_1_0.png",
                width = 1011*2,height = 692*2, res = 120*2)
p
dev.off()





fix_for_table_new <- function(model_1, model_2, new_names = NULL){
  tmp1 <- summary(model_1, verbose = FALSE)
  tmp2 <- summary(model_2, verbose = FALSE)

  tmp1$info_all <-  rbind(round(cbind(tmp1$Incidence$coefficients,1), 3),
                          round(cbind(tmp1$Duration$coefficients,0), 3))
  tmp2$info_all <-  rbind(round(cbind(tmp2$Incidence$coefficients,1), 3),
                          round(cbind(tmp2$Duration$coefficients,0), 3))
  # tmp1$info_all <- format(tmp1$info_all, nsmall = 3)
  # tmp2$info_all <- format(tmp1$info_all, nsmall = 3)
  tmp1$info_all <- tmp1$info_all[,c(1,2,5)]
  tmp2$info_all <- tmp2$info_all[,c(1,2,5)]
  tmp1 <- cbind(format(tmp1$info_all[,c(1,2)], nsmall =3),
                format(round(exp(log(2)*tmp1$info_all[,1]),3), nsmall =3),
                paste(rownames(tmp1$info_all), tmp1$info_all[,3]))
  tmp2 <- cbind(format(tmp2$info_all[,c(1,2)], nsmall =3),
                format(round(exp(log(2)*tmp2$info_all[,1]),3), nsmall =3),
                paste(rownames(tmp2$info_all), tmp2$info_all[,3]))
  rownames_new <- unique(c(tmp1[,4], tmp2[,4]))
  rownames_new_org <- rownames_new
  tmp1_new <- matrix("", nrow = length(rownames_new), ncol = 3)
  tmp2_new <- matrix("", nrow = length(rownames_new), ncol = 3)

  tmp1_new[match(tmp1[,4], rownames_new),] <- tmp1[,1:3]
  tmp2_new[match(tmp2[,4], rownames_new),] <- tmp2[,1:3]

  rownames_new <- gsub(" 0", "", rownames_new)
  rownames_new <- gsub(" 1", "", rownames_new)

  rownames(tmp1_new) <- rownames_new
  rownames(tmp2_new) <- rownames_new
  colnames(tmp1_new) <- c("$\\hat{\\alpha}$", "SE", "$2^{\\hat{\\alpha}}$")
  colnames(tmp2_new) <- colnames(tmp1_new)

  res <- cbind(tmp1_new, tmp2_new)
  if(!is.null(new_names)){
    rownames(res) <- new_names
  }
  res <- list(res= res,
              rownames_new_org =rownames_new_org)
  old_row <- rownames(res$res)
  old_col <- colnames(res$res)
  res$res <- apply(res$res,2, FUN = function(x){
    x <- as.numeric(x)
    if(TRUE%in% x<0){
      x[x<0] <- paste("$",x[x<0],"$")
    }
    return(x)
  })
  colnames(res$res) <- old_col
  rownames(res$res) <- old_row
  new_res <- res$res
  for(i in 1:ncol(res$res)){
    new_res[,i] <- sprintf("%.3f", res$res[,i])
    new_res[,i] <- gsub(pattern = "0\\.",replacement = ".",new_res[,i])
  }
  new_res <- gsub(pattern = "-",replacement = "$-$",new_res)

  ind <- new_names %in% c("Friendship Match", "Both Female", "Friendship Match",
                          "Both Female")
  new_res[,3][ind] <- NA
  new_res[,6][ind] <- NA
  new_res[new_res == "NA"] <- NA


  return(new_res)
}
fix_for_table_alt <- function(model_1, model_2, new_names = NULL){
  tmp1 <- summary(model_1, verbose = FALSE)
  tmp2 <- summary(model_2, verbose = FALSE)

  tmp1$info_all <-  rbind(round(cbind(tmp1$Incidence$coefficients,1), 3),
                          round(cbind(tmp1$Duration$coefficients,0), 3))
  tmp2$info_all <-  rbind(round(cbind(tmp2$Incidence$coefficients,1), 3),
                          round(cbind(tmp2$Duration$coefficients,0), 3))
  rownames(tmp1$info_all) <- gsub(pattern = "log", replacement = "",x = rownames(tmp1$info_all))
  rownames(tmp1$info_all) <- gsub(pattern = "identity", replacement = "",x = rownames(tmp1$info_all))

  rownames(tmp2$info_all) <- gsub(pattern = "log", replacement = "",x = rownames(tmp2$info_all))
  rownames(tmp2$info_all) <- gsub(pattern = "identity", replacement = "",x = rownames(tmp2$info_all))
  # tmp1$info_all <- format(tmp1$info_all, nsmall = 3)
  # tmp2$info_all <- format(tmp1$info_all, nsmall = 3)
  tmp1$info_all <- tmp1$info_all[,c(1,2,5)]
  tmp2$info_all <- tmp2$info_all[,c(1,2,5)]
  tmp1 <- cbind(format(tmp1$info_all[,c(1,2)], nsmall =3),
                format(round(exp(log(2)*tmp1$info_all[,1]),3), nsmall =3),
                paste(rownames(tmp1$info_all), tmp1$info_all[,3]))
  tmp2 <- cbind(format(tmp2$info_all[,c(1,2)], nsmall =3),
                format(round(exp(log(2)*tmp2$info_all[,1]),3), nsmall =3),
                paste(rownames(tmp2$info_all), tmp2$info_all[,3]))
  rownames_new <- unique(c(tmp1[,4], tmp2[,4]))
  rownames_new_org <- rownames_new
  tmp1_new <- matrix("", nrow = length(rownames_new), ncol = 3)
  tmp2_new <- matrix("", nrow = length(rownames_new), ncol = 3)

  tmp1_new[match(tmp1[,4], rownames_new),] <- tmp1[,1:3]
  tmp2_new[match(tmp2[,4], rownames_new),] <- tmp2[,1:3]

  rownames_new <- gsub(" 0", "", rownames_new)
  rownames_new <- gsub(" 1", "", rownames_new)

  rownames(tmp1_new) <- rownames_new
  rownames(tmp2_new) <- rownames_new
  colnames(tmp1_new) <- c("$\\hat{\\alpha}$", "SE", "$2^{\\hat{\\alpha}}$")
  colnames(tmp2_new) <- colnames(tmp1_new)

  res <- cbind(tmp1_new, tmp2_new)
  if(!is.null(new_names)){
    rownames(res) <- new_names
  }
  res <- list(res= res,
              rownames_new_org =rownames_new_org)
  old_row <- rownames(res$res)
  old_col <- colnames(res$res)
  res$res <- apply(res$res,2, FUN = function(x){
    x <- as.numeric(x)
    if(TRUE%in% x<0){
      x[x<0] <- paste("$",x[x<0],"$")
    }
    return(x)
  })
  colnames(res$res) <- old_col
  rownames(res$res) <- old_row
  new_res <- res$res
  for(i in 1:ncol(res$res)){
    new_res[,i] <- sprintf("%.3f", res$res[,i])
    new_res[,i] <- gsub(pattern = "0\\.",replacement = ".",new_res[,i])
  }
  new_res <- gsub(pattern = "-",replacement = "$-$",new_res)

  ind <- new_names %in% c("Friendship Match", "Both Female", "Friendship Match",
                          "Both Female")
  new_res[,3][ind] <- NA
  new_res[,6][ind] <- NA
  new_res[new_res == "NA"] <- NA


  return(new_res)
}

fix_for_table_rem <- function(model_1, model_2, new_names = NULL){
  tmp1 <- summary(model_1, verbose = FALSE)
  tmp2 <- summary(model_2, verbose = FALSE)

  tmp1$info_all <-  rbind(round(cbind(tmp1$coefficients,1), 3),
                          round(cbind(tmp1$coefficients,0), 3))
  tmp2$info_all <-  rbind(round(cbind(tmp2$coefficients,1), 3),
                          round(cbind(tmp2$coefficients,0), 3))
  # tmp1$info_all <- format(tmp1$info_all, nsmall = 3)
  # tmp2$info_all <- format(tmp1$info_all, nsmall = 3)
  tmp1$info_all <- tmp1$info_all[,c(1,2,5)]
  tmp2$info_all <- tmp2$info_all[,c(1,2,5)]
  tmp1 <- cbind(format(tmp1$info_all[,c(1,2)], nsmall =3),
                format(round(exp(log(2)*tmp1$info_all[,1]),3), nsmall =3),
                paste(rownames(tmp1$info_all), tmp1$info_all[,3]))
  tmp2 <- cbind(format(tmp2$info_all[,c(1,2)], nsmall =3),
                       format(round(exp(log(2)*tmp2$info_all[,1]),3), nsmall =3),
                paste(rownames(tmp2$info_all), tmp2$info_all[,3]))
  rownames_new <- unique(c(tmp1[,4], tmp2[,4]))
  rownames_new_org <- rownames_new
  tmp1_new <- matrix("", nrow = length(rownames_new), ncol = 3)
  tmp2_new <- matrix("", nrow = length(rownames_new), ncol = 3)

  tmp1_new[match(tmp1[,4], rownames_new),] <- tmp1[,1:3]
  tmp2_new[match(tmp2[,4], rownames_new),] <- tmp2[,1:3]

  rownames_new <- gsub(" 0", "", rownames_new)
  rownames_new <- gsub(" 1", "", rownames_new)

  rownames(tmp1_new) <- rownames_new
  rownames(tmp2_new) <- rownames_new
  colnames(tmp1_new) <- c("$\\hat{\\alpha}$", "SE", "$2^{\\hat{\\alpha}}$")
  colnames(tmp2_new) <- colnames(tmp1_new)

  res <- cbind(tmp1_new, tmp2_new)
  if(!is.null(new_names)){
    rownames(res) <- new_names
  }
  res <- list(res= res,
       rownames_new_org =rownames_new_org)
  old_row <- rownames(res$res)
  old_col <- colnames(res$res)
  res$res <- apply(res$res,2, FUN = function(x){
    x <- as.numeric(x)
    if(TRUE%in% x<0){
      x[x<0] <- paste("$",x[x<0],"$")
    }
    return(x)
  })
  colnames(res$res) <- old_col
  rownames(res$res) <- old_row
  new_res <- res$res
  for(i in 1:ncol(res$res)){
    new_res[,i] <- sprintf("%.3f", res$res[,i])
    new_res[,i] <- gsub(pattern = "0\\.",replacement = ".",new_res[,i])
  }
  # browser()
  new_res <- gsub(pattern = "-",replacement = "$-$",new_res)

  ind <- new_names %in% c("Friendship Match", "Both Female", "Friendship Match",
                          "Both Female")
  new_res[,3][ind] <- NA
  new_res[,6][ind] <- NA
  new_res[new_res == "NA"] <- NA


  return(new_res)
}

new_names <- c("Current Common Partner","General Common Partner", "Number Interaction", "Friendship Match",
               "Both Female", "Current Interaction", "Number Interaction", "Current Common Partner","General Common Partner",
               "Friendship Match","Both Female")

res <- fix_for_table_new(model_1 = models_proximity, model_2 = models, new_names = new_names)

res_sensitivity_proximity <- fix_for_table_new(model_1 = models_proximity, model_2 = models_proximity_alt, new_names = new_names)
res_sensitivity_proximity <- res_sensitivity_proximity[,-c(3,6)]

new_names_call <- c("Number Interaction", "General Common Partner", "Friendship Match",
                    "Both Female", "Current Interaction", "Number Interaction", "General Common Partner",
                    "Friendship Match","Both Female")

res_sensitivity_call <- fix_for_table_new(model_1 = models, model_2 = models_alt, new_names = new_names_call)
# Cut the columns with the exponentiated coefficients for the sensitivity analyses, as they are not informative in this case
res_sensitivity_call <- res_sensitivity_call[,-c(3,6)]

latex(object = res,n.cgroup = c(3,3), n.rgroup = c(5,6),
      rgroup = c("Incidence ($\\hat \\alpha^{0\\rightarrow 1}$)",
                 "Duration ($\\hat \\alpha^{1\\rightarrow 0}$)"),
      cgroup = c("Physical (Co-location)", "Digital (Call)"),
      file="Application/Results/table_res.tex",
      label = "tbl:results",
      title = "",
      first.hline.double= F,
      rowlabel.just = c("l"),
      caption = paste("Parameter estimates $(\\hat{\\alpha})$ and standard errors (SE) of the Durational Event Model applied to the Physical (Co-location) and call data. The column $2^{\\hat \\alpha}$ shown for $\\log(\\cdot +1)$ transformed statistics represents the effect of the first change in the respective statistic. "),
      cgroup.just="c",where = "t!")

latex(object = res_sensitivity_call,n.cgroup = c(2,2), n.rgroup = c(4,5),
      rgroup = c("Incidence ($\\hat \\alpha^{0\\rightarrow 1}$)",
                 "Duration ($\\hat \\alpha^{1\\rightarrow 0}$)"),
      cgroup = c("Every Hour", "Every two hours"),
      file="Application/Results/table_res_sensitivity_call.tex",
      label = "tbl:results_sens_call",
      title = "",
      first.hline.double= F,
      rowlabel.just = c("l"),
      caption = paste("Parameter estimates $(\\hat{\\alpha})$ and standard errors (SE) of the Durational Event Model applied to the call data assuming a baseline intensity that can change every hour (left column) and every second hour(right column). "),
      cgroup.just="c",where = "t!")

latex(object = res_sensitivity_proximity,n.cgroup = c(2,2), n.rgroup = c(5,6),
      rgroup = c("Incidence ($\\hat \\alpha^{0\\rightarrow 1}$)",
                 "Duration ($\\hat \\alpha^{1\\rightarrow 0}$)"),
      cgroup = c("Every Hour", "Every two hours"),
      file="Application/Results/table_res_sensitivity_proximity.tex",
      label = "tbl:results_sens_proximity",
      title = "",
      first.hline.double= F,
      rowlabel.just = c("l"),
      caption = paste("Parameter estimates $(\\hat{\\alpha})$ and standard errors (SE) of the Durational Event Model applied to the Physical (Co-location) data assuming a baseline intensity that can change every hour (left column) and every second hour(right column). "),
      cgroup.just="c",where = "t!")


#  Comparison of the baseline estimates

time_call <- as.numeric(gsub(pattern = "time_", "", names(models$model_0_1$est_time)))
time_prox <- as.numeric(gsub(pattern = "time_", "", names(models_proximity$model_0_1$est_time)))
length(time_call)
length(time_prox)

time_alt <- as.numeric(gsub(pattern = "time_", "", names(models_proximity_alt$model_0_1$est_time)))
time_alt <- time_alt/60/60/24 #transform time from seconds to days
time <- as.numeric(gsub(pattern = "time_", "", names(models_proximity$model_0_1$est_time)))
time <- time/60/60/24 #transform time from seconds to days

df_hourly <- data.frame(Time = c(0,time),
                        Value =  2*mean(models_proximity$model_0_1$est_degree) + c(0,models_proximity$model_0_1$est_time), Series = "Every Hour")
df_bi_hourly <- data.frame(Time = c(0,time_alt),
                           Value =  2*mean(models_proximity_alt$model_0_1$est_degree) + c(0,models_proximity_alt$model_0_1$est_time), Series = "Every Two Hours")
df_bi_hourly_0_1 <- df_bi_hourly
df_0_1 <- rbind(df_hourly, df_bi_hourly)

time_alt <- as.numeric(gsub(pattern = "time_", "", names(models_proximity_alt$model_1_0$est_time)))
time_alt <- time_alt/60/60/24 #transform time from seconds to days


time <- as.numeric(gsub(pattern = "time_", "", names(models_proximity$model_1_0$est_time)))
time <- time/60/60/24 #transform time from seconds to days

df_hourly <- data.frame(Time = c(0,time),
                        Value = 2*mean(models_proximity$model_1_0$est_degree) + c(0,models_proximity$model_1_0$est_time),
                        Series = "Every Hour")
df_bi_hourly <- data.frame(Time = c(0,time_alt),
                           Value = 2*mean(models_proximity_alt$model_1_0$est_degree) + c(0,models_proximity_alt$model_1_0$est_time), Series = "Every Two Hours")
df_1_0 <- rbind(df_hourly, df_bi_hourly)



# Plot
a <- ggplot(df_0_1, aes(x = Time, y = Value, color = Series)) +
  geom_step(size = 1) +
  geom_point(data = df_bi_hourly_0_1, aes(x = Time, y = Value), shape = 21, size = 3, fill = "white") +
  ylab(expression(f(t,hat(gamma)^{0 %->% 1}))) +
  xlab("time in days (t)") +
  my_theme +
  theme(legend.position = "top")


b <- ggplot(df_1_0, aes(x = Time, y = Value, color = Series)) +
  geom_step(size = 1) +
  geom_point(data = df_bi_hourly, aes(x = Time, y = Value), shape = 21, size = 3, fill = "white") +
  ylab(expression(f(t,hat(gamma)^{1 %->% 0}))) +
  xlab("time in days (t)") +
  my_theme+
  theme(legend.position = "none")

library(patchwork)



time_alt <- as.numeric(gsub(pattern = "time_", "", names(models_alt$model_0_1$est_time)))
time_alt <- time_alt/60/60/24 #transform time from seconds to days
time <- as.numeric(gsub(pattern = "time_", "", names(models$model_0_1$est_time)))
time <- time/60/60/24 #transform time from seconds to days

df_hourly <- data.frame(Time = c(0,time),
                        Value = 2*mean(models$model_0_1$est_degree) + c(0,models$model_0_1$est_time),
                        Series = "Every Hour")
df_bi_hourly <- data.frame(Time = c(0,time_alt),
                           Value =2*mean(models_alt$model_0_1$est_degree) + c(0,models_alt$model_0_1$est_time),
                           Series = "Every Two Hours")
df_bi_hourly_0_1 <- df_bi_hourly
df_0_1 <- rbind(df_hourly, df_bi_hourly)

time_alt <- as.numeric(gsub(pattern = "time_", "", names(models_alt$model_1_0$est_time)))
time_alt <- time_alt/60/60/24 #transform time from seconds to days


time <- as.numeric(gsub(pattern = "time_", "", names(models$model_1_0$est_time)))
time <- time/60/60/24 #transform time from seconds to days

df_hourly <- data.frame(Time = c(0,time),
                        Value = 2*mean(models$model_1_0$est_degree) +  c(0,models$model_1_0$est_time),
                        Series = "Every Hour")
df_bi_hourly <- data.frame(Time = c(0,time_alt),
                           Value = 2*mean(models_alt$model_1_0$est_degree) +  c(0,models_alt$model_1_0$est_time), Series = "Every Two Hours")
df_1_0 <- rbind(df_hourly, df_bi_hourly)

# Plot
c <- ggplot(df_0_1, aes(x = Time, y = Value, color = Series)) +
  geom_step(size = 1) +
  geom_point(data = df_bi_hourly_0_1, aes(x = Time, y = Value), shape = 21, size = 3, fill = "white") +
  ylab(expression(f(t,hat(gamma)^{0 %->% 1}))) +
  xlab("time in days (t)") +
  my_theme +
  theme(legend.position = "none")
c
d <- ggplot(df_1_0, aes(x = Time, y = Value, color = Series)) +
  geom_step(size = 1) +
  geom_point(data = df_bi_hourly, aes(x = Time, y = Value), shape = 21, size = 3, fill = "white") +
  ylab(expression(f(t,hat(gamma)^{1 %->% 0}))) +
  xlab("time in days (t)") +
  my_theme+
  theme(legend.position = "none")

add_1 <- scale_color_manual("",values = c("Every Hour" = "black", "Every Two Hours" = "grey"))
# add_2 <- scale_linetype_manual("",values = c("Every Hour" = "solid", "Every Two Hours" = "dotted"))


col_title_left <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Incidence", size = 10) +
  theme_void()

col_title_right <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Duration", size = 10) +
  theme_void()

row_title_top <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Physical (Co-location)", size = 10, angle = 90) +
  theme_void()

row_title_bottom <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Digital (Call)", size = 10, angle = 90) +
  theme_void()

a <- a + theme(plot.title = element_blank()) + add_1
b <- b + theme(plot.title = element_blank())+ add_1
c <- c + theme(plot.title = element_blank())+ add_1
d <- d + theme(plot.title = element_blank())+ add_1

design_grid <- "
  #AB
  CDE
  FGH
"

final_plot <- wrap_plots(
  A = col_title_left,   B = col_title_right,
  C = row_title_top,    D = a,               E = b,
  F = row_title_bottom, G = c,               H = d,
  design = design_grid
) +
  plot_layout(
    guides = "collect",
    widths = c(0.1, 1, 1),
    heights = c(0.1, 1, 1)
  ) &
  theme(legend.position = "top")

Cairo::CairoPNG("Application/Plots/baseline_sens.png",
                width = 1011 * 4, height = 692 * 4, res = 120 * 2)
print(final_plot)
dev.off()



# GOF Assessment ----

models <- readRDS(file = "Application/Results/gof_results_denmark.rds")
models_alt <- readRDS(file = "Application/Results/gof_results_denmark_no_endo.rds")
models_glm <- readRDS(file = "Application/Results/gof_results_denmark_no_pop.rds")

new_names <- c( "Number Interaction", "General Common Partner","Friendship Match",
               "Both Female", "Current Interaction", "Number Interaction","General Common Partner",
               "Friendship Match","Both Female")

res_1 <- fix_for_table_alt(model_1 = models, model_2 = models_alt, new_names = new_names)
res_2 <- fix_for_table_alt(model_1 = models, model_2 = models_glm, new_names = new_names)
res <- cbind(res_1[,c(1,2,4,5)], res_2[,c(4,5)])
AIC <- c(summary(models)$AIC, "",summary(models_alt)$AIC,"", summary(models_glm)$AIC,"")
res <- rbind(res, AIC)

latex(object = res,n.cgroup = c(2,2,2), n.rgroup = c(4,5,1),
      rgroup = c("Incidence ($\\hat \\alpha^{0\\rightarrow 1}$)",
                 "Duration ($\\hat \\alpha^{1\\rightarrow 0}$)",
                 "AIC"),
      cgroup = c("Full", "Reduced", "No Popularity"),
      file="Application/Results/table_res_physical.tex",
      label = "tbl:results_gof",
      title = "",
      first.hline.double= F,
      rowlabel.just = c("l"),
      caption = paste("Parameter estimates $(\\hat{\\alpha})$ and standard errors (SE) of the Durational Event Model applied to the call data."),
      cgroup.just="c",where = "t!")



end_train <- round(nrow(events_cut)*0.8)
start_test <- round(nrow(events_cut)*0.8) + 1
res_tmp <- get_ranking(models, edgelist_test = events_cut[start_test:4152,],
                       verbose = TRUE, k_max = 500)
res_glm <- get_ranking(models_glm, edgelist_test = events_cut[start_test:4152,],
                       verbose = TRUE, k_max = 500)
res_no_endo <- get_ranking(models_alt, edgelist_test = events_cut[start_test:4152,],verbose = TRUE, k_max = 500)

models_proximity <- readRDS(file = "Application/Results/gof_results_denmark_proximity.rds")
models_proximity_no_endo <- readRDS(file = "Application/Results/gof_results_denmark_no_endo_proximity.rds")
models_proximity_no_degree <- readRDS(file = "Application/Results/gof_results_denmark_no_pop_proximity.rds")

res <- summary(models_proximity)
res_alt <- summary(models_proximity_no_endo)
new_names <- c("Current Common Partner","General Common Partner", "Number Interaction", "Friendship Match",
               "Both Female", "Current Interaction", "Number Interaction", "Current Common Partner","General Common Partner",
               "Friendship Match","Both Female")
res_1 <- fix_for_table_alt(model_1 = models_proximity, model_2 = models_proximity_no_endo, new_names = new_names)

res_2 <- fix_for_table_alt(model_1 = models_proximity, model_2 = models_proximity_no_degree, new_names = new_names)
res <- cbind(res_1[,c(1,2,4,5)], res_2[,c(4,5)])
AIC <- c(round(summary(models_proximity)$AIC, digits = 3), "",round(summary(models_proximity_no_endo)$AIC, digits = 3),"",
         round(summary(models_proximity_no_degree)$AIC, digits = 3),"")
res <- rbind(res, AIC)


latex(object = res,n.cgroup = c(2,2,2), n.rgroup = c(5,6,1),
      rgroup = c("Incidence ($\\hat \\alpha^{0\\rightarrow 1}$)",
                 "Duration ($\\hat \\alpha^{1\\rightarrow 0}$)",
                 "AIC"),
      cgroup = c("Full", "Reduced", "No Popularity"),
      file="Application/Results/table_res_colocation.tex",
      label = "tbl:results_gof_proximity",
      title = "",
      first.hline.double= F,
      rowlabel.just = c("l"),
      caption = paste("Parameter estimates $(\\hat{\\alpha})$ and standard errors (SE) of the Durational Event Model applied to the co-location data."),
      cgroup.just="c",where = "t!")

end_train <- round(nrow(proximity_cut)*0.8)
start_test <- round(nrow(proximity_cut)*0.8) + 1

res_tmp_2 <- get_ranking(models_proximity,
                         edgelist_test = as.matrix(proximity_cut[start_test:nrow(proximity_cut),1:4]),
                         verbose = TRUE, k_max = 500)
res_alt_2 <- get_ranking(models_proximity_no_degree, edgelist_test = as.matrix(proximity_cut[start_test:nrow(proximity_cut),1:4]),
                         verbose = TRUE, k_max = 500)
res_no_endo_2 <- get_ranking(models_proximity_no_endo, edgelist_test = as.matrix(proximity_cut[start_test:nrow(proximity_cut),1:4]),
                             verbose = TRUE, k_max = 500)

plot_data <- data.table(res_tmp)
plot_data$Precision <- NULL
plot_data$Recall_no_endo <- res_no_endo$Recall
plot_data$Recall_no_popularity <- res_glm$Recall
plot_data <- melt.data.table(plot_data, id.vars = "Cutpoint")
levels(plot_data$variable) <- c("Full", "Reduced",
                                "No Popularity")



plot_data_2 <- data.table(res_tmp_2)
plot_data_2$Precision <- NULL
plot_data_2$Recall_no_endo <- res_no_endo_2$Recall

plot_data_2$Recall_no_popularity <- res_alt_2$Recall
plot_data_2 <- melt.data.table(plot_data_2, id.vars = "Cutpoint")
levels(plot_data_2$variable) <-c("Full", "Reduced",
                                 "No Popularity")

lt4 <- c(
  "Full"   = "solid",
  "Reduced"= "longdash",
  "No Popularity"= "dotdash"
)

plot_2 <- ggplot(data = plot_data_2) +
  geom_line(aes(x = Cutpoint, y = value, colour = variable, linetype = variable)) +
  labs(x = "Cutpoint (k)", y = "Recall(k)", color = " ",
       linetype = " ", caption = "(b) Physical (Co-location)") +
  theme_pubr(base_size = 20) +
  theme(
    plot.title = element_blank(),                  # no top title
    plot.caption.position = "plot",                # place in bottom margin
    plot.caption = element_text(
      hjust = 0.5, size = 20
    )
  ) +
  scale_linetype_manual(values = lt4) +
  scale_x_continuous(breaks =c(1,250, 500, 750, 1000),labels = label_clean_numbers) +
  scale_y_continuous(labels = label_clean_numbers) +
  scale_color_manual(values = c("Full" = "#0072B2", "Reduced" = "#E69F00",
                                "No Popularity" = "#D55E00"))




plot_2

plot_1 <- ggplot(data = plot_data) +
  geom_line(aes(x = Cutpoint, y = value, colour = variable, linetype = variable)) +
  labs(x = "Cutpoint (k)", y = "Recall(k)", linetype  = " ", color = " ", caption = "(c) Digital (Call)") +
  theme_pubr(base_size = 20) +
  scale_x_continuous(breaks =c(1,250, 500, 750, 1000),labels = label_clean_numbers) +
  scale_y_continuous(labels = label_clean_numbers) +
  theme(
    plot.title = element_blank(),                  # no top title
    plot.caption.position = "plot",                # place in bottom margin
    plot.caption = element_text(
      hjust = 0.5, size = 20
    )
  ) +
  scale_linetype_manual(values = lt4) +
  scale_color_manual(values = c("Full" = "#0072B2", "Reduced" = "#E69F00",
                                "No Popularity" = "#D55E00"))

plot_1




label_clean_numbers <- function(x) {
  sapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    if (abs(v) >= 1000) return(comma(v, accuracy = 1))
    if (v %% 1 == 0) return(formatC(v, format = "f", digits = 0))
    s <- number(v, accuracy = 0.1, big.mark = ",", trim = TRUE)
    sub("^(-?)0\\.", "\\1.", s)
  })
}

a <- ggplot(residuals, aes(x = time)) +
  geom_step(aes(y = surv, color = "Digital"), linewidth = 1) +
  geom_step(aes(y = surv, color = "Physical"), linewidth = 1, data = residuals_proximity) +
  # geom_step(aes(y = surv_proximity, color = "Empirical"), linewidth = 1) +
  geom_line(aes(y = theoretical, color = "Theoretical"),
            linewidth = 1, linetype = "dashed") +
  scale_color_manual("",values = c("Digital" = "black",
                                   "Physical" = "steelblue",
                                   "Theoretical" = "grey")) +
  theme_pubr(base_size = 20) +
  theme(
    plot.title = element_blank(),
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0.5, size = 20)
  ) +
  labs(y = TeX(r'(Survival function $P(X > t)$)'),
       x = "t",  caption = "(a) Goodness-of-Fit")+
  scale_y_continuous(labels = label_clean_numbers) +
  xlim(c(0,10))

# scale_x_discrete(labels = label_clean_numbers(unique(plot_data$n_nodes))) +


a + ((plot_2 + plot_1 + plot_layout(guides = "collect")) & theme(legend.position = "top")) +
  plot_layout(widths = c(1,2))


ggsave(filename = "Application/Plots/denmark_gof.png", width = 14, height = 5)


# Application to the REM --------------------------------------------

names(sms) <- c("timestamp", "from", "to")
sms$tmp = sms$to<sms$from
tmp = sms$to
sms$to[sms$tmp] = sms$from[sms$tmp]
sms$from[sms$tmp] = tmp[sms$tmp]
sms$tmp = NULL

actor_data$exclude <- !actor_data$id %in% unique(as.numeric(cbind(sms$from,sms$to)))
actor_data$exclude[is.na(actor_data$female)] <- TRUE
actor_data$exclude[which(rowSums(friendship_mat) ==0)] <- TRUE


sms_cut = sms[sms$from %in% actor_data$id[!actor_data$exclude] &
                sms$to %in% actor_data$id[!actor_data$exclude]]

actor_data$exclude[which(!actor_data$id %in% unique(as.numeric(cbind(sms_cut$from,sms_cut$to))))] <- TRUE
friendship_mat_cut <- friendship_mat[!actor_data$exclude, !actor_data$exclude]
both_female_cut <- both_female[!actor_data$exclude, !actor_data$exclude]
one_female_cut <- one_female[!actor_data$exclude, !actor_data$exclude]

actor_data$id_included <- NA
actor_data$id_included[!actor_data$exclude] <- 1:sum(!actor_data$exclude)

sms_cut$from <- actor_data$id_included[match(sms_cut$from,actor_data$id)]
sms_cut$to <- actor_data$id_included[match(sms_cut$to,actor_data$id)]
sort(unique(c(sms_cut$from, sms_cut$to)))

sms_cut$type <- 1
sms_cut_end <- sms_cut
sms_cut_end$type <- 0
sms_cut_end$timestamp <- sms_cut_end$timestamp + 0.1

sms_cut_all <- rbind(sms_cut, sms_cut_end)
sms_cut_all <- sms_cut_all[order(sms_cut_all$timestamp),]

models_sms <- readRDS("Application/Results/denmark_models_sms.rds")
new_names <- c("Common Partner", "Number Interaction", "Friendship Match","Both Female")
# debugonce(fix_for_table_rem)
res <- fix_for_table_rem(model_1 = models_sms, model_2 = models_sms,
                         new_names = rep(new_names,2))

latex(object = res[1:length(new_names),1:3],
      file="Application/Results/table_res_sms.tex",
      label = "tbl:results_sms",
      title = "",
      first.hline.double= F,
      rowlabel.just = c("l"),
      caption = paste("Parameter estimates $(\\hat{\\alpha})$ and standard errors (SE) of the Relational Event Model applied to the Texting data. "),
      cgroup.just="c",where = "t!")


time <- as.numeric(gsub(pattern = "time_", "", names(models_sms$model$est_time)))
time <- time/60/60/24 #transform time from seconds to days

step_data_1 <- data.frame(
  x = c(0,time),
  y = c(0,models_sms$model$est_time)
)

smoothed <- loess(step_data_1$y~ step_data_1$x)
smooth_data_2 <- data.frame(
  x = step_data_1$x,
  y = smoothed$fitted
)

step_data_2 <- data.frame(
  x = step_data_1$x,
  y = smoothed$residuals
)


a <- ggplot(step_data_1, aes(x = x, y = y)) +
  geom_step(direction = "hv", colour = "gray")  +
  geom_smooth(se = FALSE,colour = "black")+
  theme_pubr(base_size = 20) +
  ylim(c(range(step_data_1$y)))+
  ylab(expression(f(t,hat(gamma)))) +
  xlab("time in days (t)")

c <- ggplot(step_data_2, aes(x = x, y = y)) +
  geom_step(colour = "gray")  +
  # ylim(c(range(step_data_1$y)))+
  theme_pubr(base_size = 20) +
  ylab("Residual of trend component") +
  xlab("time in days (t)")

Cairo::CairoPNG("Application/Plots/baseline_texting.png",
                width = 1011*3,height = 692*2, res = 120*2)

a/c
dev.off()
n_nodes <- nrow(friendship_mat_cut)
edges <- sms_cut_all[, .(weight = log(.N)),by = .(from, to)]
edges$from <- as.character(edges$from)
edges$to <- as.character(edges$to)
# edges$weight <- log(edges$weight)
edges <- data.frame(edges)
# plot(log(models_sms$model_0_1$est_degree - min(models_sms$model_0_1$est_degree) +1))

nodes <- data.table(id = as.character(1:n_nodes), size = models_sms$model$est_degree - min(models_sms$model$est_degree))
nodes$size <- log(models_sms$model$est_degree - min(models_sms$model$est_degree) +1)*3
library(igraph)
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
set.seed(123)
Cairo::CairoPNG("Application/Plots/texting_popularity.png",
                width = 692*4,height = 692*4, res = 120*2)

plot(
  g,
  edge.width = E(g)$weight,  # Scale edge width by weight
  vertex.size = V(g)$size,
  vertex.label = NA,
  edge.color = "gray",
  vertex.color = "black"
)



dev.off()

