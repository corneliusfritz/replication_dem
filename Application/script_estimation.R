# script_estimation.R
# This script prepares network covariates, aligns node cohorts between datasets,
# specifies and runs the Durational Event Models (DEM) for calls and proximity,
# and fits the standard Relational Event Model (REM) for SMS messages.

library(data.table)
library(redeem)
library(stringr)

# 1. Load Call Logs and Node Data
# setwd to the current location of the file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
calls = fread("Application/Data/calls.csv")
actor_data = fread("Application/Data/actor_data.csv")

# Align call caller/callee IDs to the master unique actor list
calls$caller <- match(calls$caller, actor_data$unique_actors)
calls$callee <- match(calls$callee, actor_data$unique_actors)

# Keep only valid calls between study participants with duration > 0
calls <- calls[!is.na(caller) & !is.na(callee) & caller != callee & duration > 0]
names(calls) <- c("timestamp", "from", "to", "duration")

# Scale call timestamps so the first call starts at t = 0
calls$time <- calls$timestamp - min(calls$timestamp)

# Enforce undirected dyad indexing (from < to)
calls$tmp <- calls$to < calls$from
tmp <- calls$to
calls$to[calls$tmp] <- calls$from[calls$tmp]
calls$from[calls$tmp] <- tmp[calls$tmp]
calls$tmp <- NULL
calls$from_to <- paste(calls$from, calls$to, sep = "_")
calls <- calls[,.(time, duration, from, to, from_to)]


# 2. Convert Call Records to Durational Boundary Events (tmp_fun)
# This function converts raw call logs (each with a start time and duration) into
# individual start (type=1, incidence) and end (type=0, dissolution) event markers.
tmp_fun <- function(x) {
  start <- x$time
  end <- x$time + x$duration

  # Bind start and end timestamps together into a single chronological sequence
  res <- data.frame(rbind(cbind(x$from, x$to, start, 1),
                          cbind(x$from, x$to, end, 0)))
  colnames(res) <- c("from", "to", "time", "type")
  res <- res[order(res$time),]

  # Check for and resolve consecutive duplicates/overlapping intervals
  # (e.g. type difference of 0 indicates back-to-back starts or ends)
  if(length(which(diff(res$type) == 0)) > 0){
    number_errors <- length(which(diff(res$type) == 0)) / 2

    errors <- which(diff(res$type) == 0)
    exclude <- c(errors, errors+1)

    exclude <- c()
    for(i in seq(from = 1, to = length(errors), by = 2)){
      exclude <- c(exclude, errors[i]:(errors[i+1]+1))
    }

    res <- rbind(res[errors[seq(1, length(errors), 2)],],
                 res[errors[seq(2, length(errors), 2)] + 1,],
                 res[-exclude,])

    res <- res[order(res$time),]
  }
  return(res)
}

# Apply the call-to-event boundary conversion per dyad
res <- calls[, .(mat = list(tmp_fun(.SD))), by = from_to]
res <- data.table(do.call(what = rbind, res$mat))

events <- res[, .(time, from, to, type)]
events$time <- as.numeric(events$time)
events$time <- events$time - min(events$time)

colnames(actor_data) <- c("id", "id_old")


# 3. Load Covariates (Bluetooth proximity, genders, Facebook friends)
proximity <- fread("Application/Data/events.csv")
gender <- fread("Application/Data/genders.csv")
actor_data$female <- gender$female[match(actor_data$id_old, gender$`# user`)]
sms <- fread("Application/Data/sms.csv")

# Load Facebook friendship network
fb_friends <- fread("Application/Data/fb_friends.csv")
names(fb_friends) <- c("from", "to")
fb_friends <- fb_friends[!from == to]
fb_friends <- fb_friends[!duplicated(fb_friends)]

# Exclude network links outside the main participant study cohort
fb_friends$included <- (fb_friends$from %in% actor_data$id_old) & (fb_friends$to %in% actor_data$id_old)
fb_friends <- fb_friends[included == TRUE]
fb_friends$from_id <- match(fb_friends$from, actor_data$id_old)
fb_friends$to_id <- match(fb_friends$to, actor_data$id_old)

# Build the symmetric Facebook adjacency matrix
friendship_mat <- matrix(0, nrow = nrow(actor_data), ncol = nrow(actor_data))
friendship_mat[cbind(fb_friends$from_id, fb_friends$to_id)] = 1
friendship_mat[cbind(fb_friends$to_id, fb_friends$from_id)] = 1

# Gender homophily dyadic covariates
both_female <- outer(actor_data$female, actor_data$female, FUN = function(x, y) { x == y })
one_female <- outer(actor_data$female, 1 - actor_data$female, FUN = function(x, y) { x == y })


# 4. Iterative Node cohort Pruning
# We must ensure that the modeled actors are consistently present and active across
# multiple datasets (calls, bluetooth proximity, facebook friendship, gender metadata).
# We run a 10-iteration filter to discard inactive nodes or nodes with missing covariate values.
actor_data$exclude <- !actor_data$id %in% unique(as.numeric(cbind(events$from, events$to)))
actor_data$exclude[is.na(actor_data$female)] <- TRUE
actor_data$exclude[which(rowSums(friendship_mat) == 0)] <- TRUE

for(i in 1:10){
  events_cut <- events[events$from %in% actor_data$id[!actor_data$exclude] &
                       events$to %in% actor_data$id[!actor_data$exclude]]

  actor_data$exclude[which(!actor_data$id %in% unique(as.numeric(cbind(events_cut$from, events_cut$to))))] <- TRUE

  proximity_cut <- proximity[proximity$from %in% actor_data$id[!actor_data$exclude] &
                             proximity$to %in% actor_data$id[!actor_data$exclude]]

  actor_data$exclude[which(!actor_data$id %in% unique(as.numeric(cbind(proximity_cut$from, proximity_cut$to))))] <- TRUE
  cat(sum(!actor_data$exclude), "\n")
}

# Subset covariate matrices for the final active cohort
friendship_mat_cut <- friendship_mat[!actor_data$exclude, !actor_data$exclude]
both_female_cut <- both_female[!actor_data$exclude, !actor_data$exclude]
one_female_cut <- one_female[!actor_data$exclude, !actor_data$exclude]

# Re-map actor IDs to a continuous sequence from 1 to N_nodes for the active cohort
actor_data$id_included <- NA
actor_data$id_included[!actor_data$exclude] <- 1:sum(!actor_data$exclude)
events_cut$from <- actor_data$id_included[match(events_cut$from, actor_data$id)]
events_cut$to <- actor_data$id_included[match(events_cut$to, actor_data$id)]
sort(unique(c(events_cut$from, events_cut$to)))

# 5. Define Time Changepoints for Non-parametric Baseline Intensity
# Set hourly changepoints (3600 seconds) for baseline baseline(changepoints)
time_changepoints <- seq(from = 60*60*24, by = 60*60*24, to = max(proximity$time))
time_changepoints <- time_changepoints[-length(time_changepoints)]
labels_changepoints <- c(1:(length(time_changepoints)))
events_cut <- events_cut[order(time)]
events_cut <- as.matrix(events_cut)
n_nodes <- nrow(friendship_mat_cut)

cp_seq <- seq(from = 3600, to = max(events_cut[, 1]), by = 3600)


# 6. Specify Durational Event Model (DEM) Formulas
formula_0_1 <-  ~
  degree +
  inertia(transformation = "log") +
  general_common_partners(transformation = "log") +
  dyadic_cov(data = friendship_mat_cut) +
  dyadic_cov(data = both_female_cut) +
  baseline(changepoints = cp_seq)

# formula_1_0: Dissolution rate (ending an active interaction)
formula_1_0 <-  ~
  degree +
  current_interaction(transformation = "log") +
  inertia(transformation = "log") +
  general_common_partners(transformation = "log") +
  dyadic_cov(data = friendship_mat_cut) +
  dyadic_cov(data = both_female_cut) +
  baseline(changepoints = cp_seq)

options(expressions = 500000)


# 7. Model Estimation: Digital Interaction (Calls)
# Fit the Durational Event Model (DEM) with hourly baseline changepoints
models <- dem(events = events_cut,
              formula_1_0 = formula_1_0,
              formula_0_1 = formula_0_1,
              n_nodes = n_nodes,
              directed = FALSE,
              simultaneous_interactions = FALSE,
              semiparametric = FALSE,estimate_0_1 = F,
              control = control.redeem(accelerated = FALSE,
                                       verbose = TRUE, return_data = TRUE,
                                       weighting = FALSE,
                                       it_max = c(2000, 10000), tol = 0.001))

# Sensitivity Analysis: Fit model with 2-hour changepoints (cp_seq_alt)
cp_seq_alt <- seq(from = 2*3600, to = max(events_cut[,1]), by = 2*3600)
formula_0_1_alt <-  ~ degree +
  inertia(transformation = "log") +
  general_common_partners(transformation = "log") +
  dyadic_cov(data = friendship_mat_cut)  +  dyadic_cov(data = both_female_cut) +
  baseline(changepoints = cp_seq)

formula_1_0_alt <-  ~ degree + current_interaction(transformation = "log") +
  inertia(transformation = "log") +
  general_common_partners(transformation = "log") +
  dyadic_cov(data = friendship_mat_cut)  +  dyadic_cov(data = both_female_cut) +
  baseline(changepoints = cp_seq)

models_alt <- dem(events = events_cut,
                  formula_1_0 = formula_1_0_alt,
                  formula_0_1 = formula_0_1_alt,
                  simultaneous_interactions = FALSE,
                  n_nodes = n_nodes, semiparametric = FALSE,
                  control = control.redeem(accelerated = FALSE,
                                           verbose = TRUE,
                                           weighting = FALSE,
                                           it_max = c(2000, 10000), tol = 0.001))


saveRDS(models, file = "Application/Results/denmark_models_call.rds")
saveRDS(models_alt, file = "Application/Results/denmark_models_call_alt.rds")
models_alt <- readRDS(file = "Application/Results/denmark_models_call.rds")
summary(models_alt$model_1_0)
summary(models$model_1_0)
# 8. Model Estimation: Physical Interaction (Co-location)
# Apply the DEM to the continuous Bluetooth co-location data.
formula_0_1 <-  ~ degree + current_common_partners(transformation = "log") +
  general_common_partners(transformation = "log") +
  inertia(transformation = "log")  +
  dyadic_cov(data = friendship_mat_cut)  +  dyadic_cov(data = both_female_cut) +
  baseline(changepoints = cp_seq)

formula_1_0 <-  ~ degree +  current_interaction(transformation = "log")  +
  inertia(transformation = "log")  +
  current_common_partners(transformation = "log") +
  general_common_partners(transformation = "log") +
  dyadic_cov(data = friendship_mat_cut)  +  dyadic_cov(data = both_female_cut) +
  baseline(changepoints = cp_seq)

training_start <- 0
proximity_cut$from <- actor_data$id_included[match(proximity_cut$from, actor_data$id)]
proximity_cut$to <- actor_data$id_included[match(proximity_cut$to, actor_data$id)]
proximity_cut <- as.matrix(proximity_cut)

# Fit Bluetooth proximity DEM.
# simultaneous_interactions = TRUE, because multiple separate pairs
# can be in proximity at the same time.
# computing on large datasets.
models_proximity <- dem(events = as.matrix(proximity_cut[,1:4]),
                        formula_1_0 = formula_1_0,
                        formula_0_1 = formula_0_1,
                        n_nodes = n_nodes,
                        training_start = training_start,
                        semiparametric = FALSE,
                        simultaneous_interactions = TRUE,
                        control = control.redeem(accelerated = FALSE,
                                                 verbose = TRUE, return_data = TRUE,
                                                 weighting = FALSE,
                                                 subsample = 0.2,
                                                 it_max = c(2000, 4000), tol = 0.0001))

# Fit Proximity Sensitivity Model with 2-hour baseline intervals
cp_seq_alt_prox <- seq(from = 2*3600, to = max(proximity_cut[,1]), by = 2*3600)
formula_0_1_alt_prox <-  ~ degree + current_common_partners(transformation = "log") +
  general_common_partners(transformation = "log") +
  inertia(transformation = "log")  +
  dyadic_cov(data = friendship_mat_cut)  +  dyadic_cov(data = both_female_cut) +
  baseline(changepoints = cp_seq_alt_prox)

formula_1_0_alt_prox <-  ~ degree +  current_interaction(transformation = "log")  +
  inertia(transformation = "log")  +
  current_common_partners(transformation = "log") +
  general_common_partners(transformation = "log") +
  dyadic_cov(data = friendship_mat_cut)  +  dyadic_cov(data = both_female_cut) +
  baseline(changepoints = cp_seq_alt_prox)

models_proximity_alt <- dem(events = as.matrix(proximity_cut[,1:4]),
                             formula_1_0 = formula_1_0_alt_prox,
                             formula_0_1 = formula_0_1_alt_prox,
                             n_nodes = n_nodes,
                             training_start = training_start,
                             semiparametric = FALSE,
                            simultaneous_interactions = TRUE,
                             control = control.redeem(accelerated = FALSE,
                                                      verbose = TRUE,
                                                      weighting = FALSE,
                                                      subsample = 0.2,
                                                      it_max = c(2000, 4000), tol = 0.0001))

saveRDS(models_proximity, file = "Application/Results/denmark_models_proximity.rds")
saveRDS(models_proximity_alt, file = "Application/Results/denmark_models_proximity_alt.rds")


# 9. Model Estimation: Relational Event Model (REM) on Instantaneous SMS Messages
# SMS messages are instantaneous point events (no duration). We model them using the
# classic Relational Event Model (REM) via the rem() function in redeem.
names(sms) <- c("timestamp", "from", "to")
sms$tmp <- sms$to < sms$from
tmp <- sms$to
sms$to[sms$tmp] <- sms$from[sms$tmp]
sms$from[sms$tmp] <- tmp[sms$tmp]
sms$tmp <- NULL

actor_data$exclude <- !actor_data$id %in% unique(as.numeric(cbind(sms$from, sms$to)))
actor_data$exclude[is.na(actor_data$female)] <- TRUE
actor_data$exclude[which(rowSums(friendship_mat) == 0)] <- TRUE

sms_cut <- sms[sms$from %in% actor_data$id[!actor_data$exclude] &
               sms$to %in% actor_data$id[!actor_data$exclude]]

actor_data$exclude[which(!actor_data$id %in% unique(as.numeric(cbind(sms_cut$from, sms_cut$to))))] <- TRUE
friendship_mat_cut <- friendship_mat[!actor_data$exclude, !actor_data$exclude]
both_female_cut <- both_female[!actor_data$exclude, !actor_data$exclude]
one_female_cut <- one_female[!actor_data$exclude, !actor_data$exclude]

actor_data$id_included <- NA
actor_data$id_included[!actor_data$exclude] <- 1:sum(!actor_data$exclude)

sms_cut$from <- actor_data$id_included[match(sms_cut$from, actor_data$id)]
sms_cut$to <- actor_data$id_included[match(sms_cut$to, actor_data$id)]
sort(unique(c(sms_cut$from, sms_cut$to)))

time_changepoints <- unique(sms_cut$time)
sms_cut$type <- 1 # SMS events only represent incidence (type=1)

check_matrix(as.matrix(sms_cut))

formula_0_1 <-  ~ degree + general_common_partners(transformation = "log") +
  inertia(transformation = "log")  +
  dyadic_cov(data = friendship_mat_cut)  +  dyadic_cov(data = both_female_cut) +
  baseline(changepoints = seq(from = 3600, to = max(events_cut[, 1]), by = 3600))

n_nodes <- nrow(one_female_cut)
training_start <- 0
time_changepoints <- time_changepoints[-length(time_changepoints)]
labels_changepoints <- c(1:(length(time_changepoints)))

# Fit the REM model using rem()
models <- rem(events = as.matrix(sms_cut),
              formula = formula_0_1,
              n_nodes = n_nodes,
              training_start = training_start,
              control = control.redeem(accelerated = TRUE,
                                       verbose = TRUE,
                                       subsample = 1,
                                       weighting = FALSE,
                                       it_max = 2000, tol = 0.0001))

saveRDS(models, file = "Application/Results/denmark_models_sms.rds")

