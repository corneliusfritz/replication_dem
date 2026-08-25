# script_oos.R
# This script performs out-of-sample (OOS) validation and goodness-of-fit (GoF) analysis.
# It splits both call and proximity datasets into training (first 80% of events)
# and testing (remaining 20%) sets, fits three competing model specifications to the
# training data, and saves these models for OOS likelihood comparisons.

library(data.table)
library(redeem)
library(stringr)
library(ggpubr)

# 1. Load Call Logs and Node Data
calls <- fread("Application/Data/calls.csv")
actor_data <- fread("Application/Data/actor_data.csv")

# Align caller/callee IDs to study participants
calls$caller <- match(calls$caller, actor_data$unique_actors)
calls$callee <- match(calls$callee, actor_data$unique_actors)

# Filter valid calls
calls <- calls[!is.na(caller) & !is.na(callee) & caller != callee & duration > 0]
names(calls) <- c("timestamp", "from", "to", "duration")
calls$time <- calls$timestamp - min(calls$timestamp)

# Enforce undirected dyad property (from < to)
calls$tmp <- calls$to < calls$from
tmp <- calls$to
calls$to[calls$tmp] <- calls$from[calls$tmp]
calls$from[calls$tmp] <- tmp[calls$tmp]
calls$tmp <- NULL
calls$from_to <- paste(calls$from, calls$to, sep = "_")
calls <- calls[, .(time, duration, from, to, from_to)]


# 2. Convert Call Records to Event Boundaries (tmp_fun)
tmp_fun <- function(x) {
  start <- x$time
  end <- x$time + x$duration
  res <- data.frame(rbind(
    cbind(x$from, x$to, start, 1),
    cbind(x$from, x$to, end, 0)
  ))
  colnames(res) <- c("from", "to", "time", "type")
  res <- res[order(res$time), ]

  # Resolve overlaps or back-to-back state changes
  if (length(which(diff(res$type) == 0)) > 0) {
    number_errors <- length(which(diff(res$type) == 0)) / 2

    errors <- which(diff(res$type) == 0)
    exclude <- c(errors, errors + 1)

    exclude <- c()
    for (i in seq(from = 1, to = length(errors), by = 2)) {
      exclude <- c(exclude, errors[i]:(errors[i + 1] + 1))
    }

    res <- rbind(
      res[errors[seq(1, length(errors), 2)], ],
      res[errors[seq(2, length(errors), 2)] + 1, ],
      res[-exclude, ]
    )

    res <- res[order(res$time), ]
  }
  return(res)
}

res <- calls[, .(mat = list(tmp_fun(.SD))), by = from_to]
res <- data.table(do.call(what = rbind, res$mat))

events <- res[, .(time, from, to, type)]
events$time <- as.numeric(events$time)
events$time <- events$time - min(events$time)

colnames(actor_data) <- c("id", "id_old")


# 3. Load Covariates (genders, Facebook friendship network)
proximity <- fread("Application/Data/events.csv")
gender <- fread("Application/Data/genders.csv")
actor_data$female <- gender$female[match(actor_data$id_old, gender$`# user`)]
sms <- fread("Application/Data/sms.csv")

fb_friends <- fread("Application/Data/fb_friends.csv")
names(fb_friends) <- c("from", "to")
fb_friends <- fb_friends[!from == to]
fb_friends <- fb_friends[!duplicated(fb_friends)]
fb_friends$included <- (fb_friends$from %in% actor_data$id_old) & (fb_friends$to %in% actor_data$id_old)
fb_friends <- fb_friends[included == TRUE]
fb_friends$from_id <- match(fb_friends$from, actor_data$id_old)
fb_friends$to_id <- match(fb_friends$to, actor_data$id_old)

friendship_mat <- matrix(0, nrow = nrow(actor_data), ncol = nrow(actor_data))
friendship_mat[cbind(fb_friends$from_id, fb_friends$to_id)] <- 1
friendship_mat[cbind(fb_friends$to_id, fb_friends$from_id)] <- 1

both_female <- outer(actor_data$female, actor_data$female, FUN = function(x, y) {
  x == y
})
one_female <- outer(actor_data$female, 1 - actor_data$female, FUN = function(x, y) {
  x == y
})


# 4. Iterative Active Cohort Node Pruning (10-iterations)
actor_data$exclude <- !actor_data$id %in% unique(as.numeric(cbind(events$from, events$to)))
actor_data$exclude[is.na(actor_data$female)] <- TRUE
actor_data$exclude[which(rowSums(friendship_mat) == 0)] <- TRUE

for (i in 1:10) {
  events_cut <- events[events$from %in% actor_data$id[!actor_data$exclude] &
    events$to %in% actor_data$id[!actor_data$exclude]]

  actor_data$exclude[which(!actor_data$id %in% unique(as.numeric(cbind(events_cut$from, events_cut$to))))] <- TRUE

  proximity_cut <- proximity[proximity$from %in% actor_data$id[!actor_data$exclude] &
    proximity$to %in% actor_data$id[!actor_data$exclude]]

  actor_data$exclude[which(!actor_data$id %in% unique(as.numeric(cbind(proximity_cut$from, proximity_cut$to))))] <- TRUE
  cat(sum(!actor_data$exclude), "\n")
}

friendship_mat_cut <- friendship_mat[!actor_data$exclude, !actor_data$exclude]
both_female_cut <- both_female[!actor_data$exclude, !actor_data$exclude]
one_female_cut <- one_female[!actor_data$exclude, !actor_data$exclude]

actor_data$id_included <- NA
actor_data$id_included[!actor_data$exclude] <- 1:sum(!actor_data$exclude)
events_cut$from <- actor_data$id_included[match(events_cut$from, actor_data$id)]
events_cut$to <- actor_data$id_included[match(events_cut$to, actor_data$id)]


# 5. Out-of-Sample Partitioning (80% Train, 20% Test)
# We set up training and testing indices based on chronological call event sequence.
time_changepoints <- seq(from = 60 * 60 * 24, by = 60 * 60 * 24, to = max(proximity$time))
time_changepoints <- time_changepoints[-length(time_changepoints)]
labels_changepoints <- c(1:(length(time_changepoints)))
events_cut <- events_cut[order(time)]
events_cut <- as.matrix(events_cut)
n_nodes <- nrow(friendship_mat_cut)

# Define boundaries for 80% train and 20% test subsets
end_train <- round(nrow(events_cut) * 0.8)
start_test <- round(nrow(events_cut) * 0.8) + 1


# 6. Specify Goodness-of-Fit (GoF) Model Formula Variants
# We fit three models: Full model, No-Endogenous model, and No-Popularity (nd) model.

# Model A: Full DEM Model Formulas
formula_0_1 <- ~
  degree + inertia(transformation = "log") +
    general_common_partners(transformation = "log") +
    dyadic_cov(data = friendship_mat_cut) + dyadic_cov(data = both_female_cut) +
    baseline(changepoints = seq(from = 3600, to = max(events_cut[1:end_train, 1]), by = 3600))

formula_1_0 <- ~ current_interaction(transformation = "log") +
  degree + inertia(transformation = "log") +
  general_common_partners(transformation = "log") +
  dyadic_cov(data = friendship_mat_cut) + dyadic_cov(data = both_female_cut) +
  baseline(changepoints = seq(from = 3600, to = max(events_cut[1:end_train, 1]), by = 3600))

# Model B: No-Popularity/No-Degree (nd) DEM Model Formulas
formula_0_1_nd <- ~
  inertia(transformation = "identity") +
    general_common_partners(transformation = "identity") +
    dyadic_cov(data = friendship_mat_cut) + dyadic_cov(data = both_female_cut) +
    baseline(changepoints = seq(from = 3600, to = max(events_cut[1:end_train, 1]), by = 3600))

formula_1_0_nd <- ~ current_interaction(transformation = "identity") +
  inertia(transformation = "identity") +
  general_common_partners(transformation = "identity") +
  dyadic_cov(data = friendship_mat_cut) + dyadic_cov(data = both_female_cut) +
  baseline(changepoints = seq(from = 3600, to = max(events_cut[1:end_train, 1]), by = 3600))

options(expressions = 500000)


# 7. Model Fitting: Digital Interaction Training Sets (Calls)
# Fit the Full DEM model to training calls
models <- dem(
  events = events_cut[1:end_train, ],
  formula_1_0 = formula_1_0,
  formula_0_1 = formula_0_1,
  n_nodes = n_nodes,
  semiparametric = FALSE,
  simultaneous_interactions = FALSE,
  control = control.redeem(
    accelerated = FALSE,
    verbose = TRUE,
    weighting = FALSE,
    it_max = c(2000, 4000), tol = 0.001
  )
)

# Fit Model C: No-Endogenous model (removes inertia and common partners, retains dyadic network covs)
models_alt <- dem(
  events = events_cut[1:end_train, ],
  formula_1_0 = ~ degree + dyadic_cov(data = friendship_mat_cut) +
    dyadic_cov(data = both_female_cut) + baseline(changepoints = seq(from = 3600, to = max(events_cut[1:end_train, 1]), by = 3600)),
  formula_0_1 = ~ degree + dyadic_cov(data = friendship_mat_cut) + dyadic_cov(data = both_female_cut) +
    baseline(changepoints = seq(from = 3600, to = max(events_cut[1:end_train, 1]), by = 3600)),
  n_nodes = n_nodes, semiparametric = FALSE,
  simultaneous_interactions = FALSE,
  control = control.redeem(
    accelerated = FALSE,
    verbose = TRUE,
    weighting = FALSE,
    it_max = c(2000, 4000), tol = 0.001
  )
)

# Fit Model B: No-Popularity model (retains endogenous dynamics, excludes degree popularity, weighting = TRUE)
models_glm <- dem(
  events = events_cut[1:end_train, ],
  formula_1_0 = formula_1_0_nd,
  formula_0_1 = formula_0_1_nd,
  n_nodes = n_nodes, semiparametric = FALSE,
  simultaneous_interactions = FALSE,
  control = control.redeem(
    accelerated = FALSE,
    verbose = TRUE,
    weighting = TRUE,
    it_max = c(2000, 4000), tol = 0.001
  )
)

# Save training results for out-of-sample validation comparison
saveRDS(models, file = "Application/Results/gof_results_denmark.rds")
saveRDS(models_alt, file = "Application/Results/gof_results_denmark_no_endo.rds")
saveRDS(models_glm, file = "Application/Results/gof_results_denmark_no_pop.rds")


# 8. Model Fitting: Physical Interaction Training Sets (Co-location Proximity)
# Align and partition proximity co-location data.
training_start <- 0
proximity_cut$from <- actor_data$id_included[match(proximity_cut$from, actor_data$id)]
proximity_cut$to <- actor_data$id_included[match(proximity_cut$to, actor_data$id)]
proximity_cut <- as.matrix(proximity_cut)
end_train <- round(nrow(proximity_cut) * 0.8)
start_test <- round(nrow(proximity_cut) * 0.8) + 1

# Define model formulas for physical co-location training
formula_0_1 <- ~ degree + current_common_partners(transformation = "log") +
  general_common_partners(transformation = "log") +
  inertia(transformation = "log") +
  dyadic_cov(data = friendship_mat_cut) + dyadic_cov(data = both_female_cut) +
  baseline(changepoints = seq(from = 3600, to = max(proximity_cut[1:end_train, 1]), by = 3600))

formula_1_0 <- ~ degree + current_interaction(transformation = "log") +
  inertia(transformation = "log") +
  current_common_partners(transformation = "log") +
  general_common_partners(transformation = "log") +
  dyadic_cov(data = friendship_mat_cut) + dyadic_cov(data = both_female_cut) +
  baseline(changepoints = seq(from = 3600, to = max(proximity_cut[1:end_train, 1]), by = 3600))

# Proximity No-Popularity (nd) model formulas
formula_0_1_nd <- ~ current_common_partners(transformation = "identity") +
  general_common_partners(transformation = "identity") +
  inertia(transformation = "identity") +
  dyadic_cov(data = friendship_mat_cut) + dyadic_cov(data = both_female_cut) +
  baseline(changepoints = seq(from = 3600, to = max(proximity_cut[1:end_train, 1]), by = 3600))

formula_1_0_nd <- ~ current_interaction(transformation = "identity") +
  inertia(transformation = "identity") +
  current_common_partners(transformation = "identity") +
  general_common_partners(transformation = "identity") +
  dyadic_cov(data = friendship_mat_cut) + dyadic_cov(data = both_female_cut) +
  baseline(changepoints = seq(from = 3600, to = max(proximity_cut[1:end_train, 1]), by = 3600))

# Proximity No-Endogenous model formulas
formula_no_endo <- ~ degree + dyadic_cov(data = friendship_mat_cut) + dyadic_cov(data = both_female_cut) +
  baseline(changepoints = seq(from = 3600, to = max(proximity_cut[1:end_train, 1]), by = 3600))

# Fit Proximity Model A (Full model) to training co-locations
models_proximity <- dem(
  events = as.matrix(proximity_cut[1:end_train, 1:4]),
  formula_1_0 = formula_1_0,
  formula_0_1 = formula_0_1,
  n_nodes = n_nodes,
  training_start = training_start,
  semiparametric = FALSE,
  simultaneous_interactions = TRUE,
  control = control.redeem(
    accelerated = FALSE,
    verbose = TRUE,
    weighting = FALSE,
    subsample = 0.2,
    it_max = c(2000, 4000), tol = 0.0001
  )
)

# Fit Proximity Model C (No-Endogenous model)
models_proximity_no_endo <- dem(
  events = as.matrix(proximity_cut[1:end_train, 1:4]),
  formula_1_0 = formula_no_endo,
  formula_0_1 = formula_no_endo,
  n_nodes = n_nodes,
  training_start = training_start,
  semiparametric = FALSE,
  simultaneous_interactions = TRUE,
  control = control.redeem(
    accelerated = FALSE,
    verbose = TRUE,
    weighting = FALSE,
    subsample = 1,
    it_max = c(2000, 4000), tol = 0.0001
  )
)

# Fit Proximity Model B (No-Popularity model)
models_proximity_no_degree <- dem(
  events = as.matrix(proximity_cut[1:end_train, 1:4]),
  formula_1_0 = formula_1_0_nd,
  formula_0_1 = formula_0_1_nd,
  n_nodes = n_nodes,
  training_start = training_start,
  semiparametric = FALSE,
  simultaneous_interactions = TRUE,
  control = control.redeem(
    accelerated = FALSE,
    verbose = TRUE,
    weighting = FALSE,
    subsample = 0.2,
    it_max = c(2000, 4000), tol = 0.0001
  )
)


saveRDS(models_proximity, file = "Application/Results/gof_results_denmark_proximity.rds")
saveRDS(models_proximity_no_endo, file = "Application/Results/gof_results_denmark_no_endo_proximity.rds")
saveRDS(models_proximity_no_degree, file = "Application/Results/gof_results_denmark_no_pop_proximity.rds")
