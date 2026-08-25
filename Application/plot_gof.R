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

# GOF Assessment ----

models <- readRDS(file = "Application/Results/gof_results_denmark.rds")
models_alt <- readRDS(file = "Application/Results/gof_results_denmark_no_endo.rds")
models_glm <- readRDS(file = "Application/Results/gof_results_denmark_no_pop.rds")

# tmp_models <- predict(models,time = max(events_cut[1:end_train,1])-1,
#                       type = "terms")
#
# tmp_models_alt <- predict(models_glm,time = max(events_cut[1:end_train,1])-1,
#                       type = "terms")
# comparable <- which(is.finite(tmp_models$formation_degree_effects) & is.finite(tmp_models$dissolution_degree_effects))
# colSums(tmp_models[comparable,4:ncol(tmp_models)])
# colSums(tmp_models_alt[comparable,4:ncol(tmp_models_alt)])

# plot(models_glm$model_0_1$data$prediction)

end_train <- round(nrow(events_cut)*0.8)
start_test <- round(nrow(events_cut)*0.8) + 1

end_test <- nrow(events_cut)
edgelist_train = events_cut[1:start_test,]
active_nodes <- unique(as.vector(edgelist_train[,2:3]))
edgelist_test = events_cut[(start_test+1):end_test,]
# edgelist_test <- edgelist_test[edgelist_test[,2] %in% active_nodes & edgelist_test[,3] %in% active_nodes,]

# active_train <- as.numeric(edgelist_train[,c(2,3)])
# active_train <- factor(active_train, levels = 1:nrow(friendship_mat_cut))
# active_test <- as.numeric(edgelist_test[,c(2,3)])
# active_test <- factor(active_test, levels = 1:nrow(friendship_mat_cut))
# plot(as.numeric(table(active_test)), as.numeric(table(active_train)))

oos_models <- get_oos_likelihood(models, edgelist_test = edgelist_test)
oos_models_glm <- get_oos_likelihood(models_glm, edgelist_test = edgelist_test)
oos_res_no_endo <- get_oos_likelihood(models_alt, edgelist_test = edgelist_test)

# mean(oos_models[is.finite(oos_models)])
# mean(oos_models_glm[is.finite(oos_models)])

res_tmp <- get_ranking(models,
                       edgelist_train =edgelist_train,
                       edgelist_test = edgelist_test,baseline_method = "last",
                       verbose = TRUE, k_max = 500)


res_glm <- get_ranking(models_glm,
                       edgelist_train =edgelist_train,
                       edgelist_test = edgelist_test,baseline_method = "last",
                       verbose = TRUE, k_max = 500)


res_no_endo <- get_ranking(models_alt,
                           edgelist_train =edgelist_train,
                           edgelist_test = edgelist_test,baseline_method = "last",
                           verbose = TRUE, k_max = 500)


plot(res_tmp$Cutpoint,res_tmp$Recall, type = "l")
lines(res_glm$Cutpoint,res_glm$Recall,col = "red")
lines(res_no_endo$Cutpoint,res_no_endo$Recall,col = "blue")


models_proximity <- readRDS(file = "Application/Results/gof_results_denmark_proximity.rds")
models_proximity_no_endo <- readRDS(file = "Application/Results/gof_results_denmark_no_endo_proximity.rds")
models_proximity_no_degree <- readRDS(file = "Application/Results/gof_results_denmark_no_pop_proximity.rds")



end_train <- round(nrow(proximity_cut)*0.8)
start_test <- round(nrow(proximity_cut)*0.8) + 1

edgelist_train = proximity_cut[1:start_test,]
edgelist_train <- as.matrix(edgelist_train)
active_nodes <- unique(as.vector(edgelist_train[,2:3]))
edgelist_test = proximity_cut[(start_test+1):nrow(proximity_cut),]
edgelist_test <- as.matrix(edgelist_test)
edgelist_test <- edgelist_test[edgelist_test[,2] %in% active_nodes & edgelist_test[,3] %in% active_nodes,]

active_train <- as.numeric(edgelist_train[,c(2,3)])
active_train <- factor(active_train, levels = 1:nrow(friendship_mat_cut))
active_test <- as.numeric(edgelist_test[,c(2,3)])
active_test <- factor(active_test, levels = 1:nrow(friendship_mat_cut))
plot(as.numeric(table(active_test)), as.numeric(table(active_train)))


# debugonce(get_ranking)
# models_proximity$events <- proximity_cut[1,]
# models_proximity$model_0_1$est_degree <- rep(mean(models_proximity$model_0_1$est_degree), length(models_proximity$model_0_1$est_degree))
# models_proximity$model_1_0$est_degree <- rep(mean(models_proximity$model_1_0$est_degree), length(models_proximity$model_1_0$est_degree))
# debugonce(get_ranking)
res_tmp_2 <- get_ranking(models_proximity,
                         edgelist_train =edgelist_train,
                         edgelist_test = edgelist_test,
                         verbose = TRUE, k_max = 500, baseline_method = "last")


res_alt_2 <- get_ranking(models_proximity_no_degree,
                         edgelist_train =edgelist_train,
                         edgelist_test = edgelist_test,
                         verbose = TRUE, k_max = 500, baseline_method = "last")

res_no_endo_2 <- get_ranking(models_proximity_no_endo,
                             edgelist_train =edgelist_train,
                             edgelist_test = edgelist_test,
                             verbose = TRUE, k_max = 500, baseline_method = "last")

# models_proximity_no_endo$model_0_1$est_time[length(models_proximity_no_endo$model_0_1$est_time)]
# models_proximity_no_endo$model_1_0$est_time[length(models_proximity_no_endo$model_0_1$est_time)]

plot(res_tmp_2$Cutpoint,res_tmp_2$Recall, type = "l")
lines(res_alt_2$Cutpoint,res_alt_2$Recall,col = "red")
lines(res_no_endo_2$Cutpoint,res_no_endo_2$Recall,col = "blue")

