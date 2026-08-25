run_simulation_for_n_node_time <- function(seed, n_nodes, simultaneous_interactions, it_max, accelerated, estimate, K, verbose, clust, path) {
  set.seed(seed)
  continuous_cov <- rnorm(n = n_nodes)
  continuous_cov <- outer(X = continuous_cov, Y = continuous_cov, FUN = function(x, y) {
    abs(x - y)
  })
  categorical_cov <- sample(x = 1:3, size = n_nodes, replace = T)
  categorical_cov <- outer(X = categorical_cov, Y = categorical_cov, FUN = function(x, y) {
    x == y
  })

  popularity_0_1 <- rnorm(n = n_nodes, sd = 1) - (6 + 0.1 * log(n_nodes))
  popularity_1_0 <- rnorm(n = n_nodes, sd = 1) - (1.6 + 0.1 * log(n_nodes))


  time_changepoints <- seq(0, 10000, length.out = 10)[-c(1, 10)]

  baseline_0_1_gt <- -seq(from = 0, to = -0.1, length.out = length(time_changepoints) + 1)[-1]
  baseline_1_0_gt <- seq(from = 0, to = -0.1, length.out = length(time_changepoints) + 1)[-1]
  formula_0_1 <- ~ degree + current_common_partners(transformation = "log") +
    dyadic_cov(data = continuous_cov) + dyadic_cov(data = categorical_cov)
  coef_0_1 <- c(
    "current_common_partners" = -0.5,
    "continuous_cov" = 1,
    "categorical_cov" = 0.5
  )
  formula_1_0 <- ~ degree + inertia(transformation = "log") +
    dyadic_cov(data = continuous_cov) + dyadic_cov(data = categorical_cov)
  coef_1_0 <- c(
    "inertia" = 0.5,
    "continuous_cov" = 0.5, "categorical_cov" = 0.5
  )
  # Load the continuous and categorical variable from the working directory to the cluster
  clusterExport(clust,
    varlist = c("continuous_cov", "categorical_cov"),
    envir = environment()
  )

  # debugonce(simulation_compare_time_memory)
  # simulation_compare_time_memory(1,formula_0_1 = formula_0_1, formula_1_0 = formula_1_0,
  #                                coef_0_1 = coef_0_1, coef_1_0 = coef_1_0,
  #                              #                                baseline_0_1_gt = baseline_0_1_gt,
  #                                baseline_1_0_gt = baseline_1_0_gt,
  #                                n_nodes = n_nodes, popularity_0_1 = popularity_0_1,
  #                                popularity_1_0 = popularity_1_0,estimate=estimate,
  #                                simultaneous_interactions = simultaneous_interactions,
  #                                verbose = FALSE, it_max = it_max,accelerated = accelerated)

  res_simulation <- parLapply(
    cl = clust, X = 1:(K), fun = simulation_compare_time_memory,
    formula_0_1 = formula_0_1, formula_1_0 = formula_1_0,
    coef_0_1 = coef_0_1, coef_1_0 = coef_1_0,
    baseline_0_1_gt = baseline_0_1_gt,
    baseline_1_0_gt = baseline_1_0_gt,
    time_changepoints = time_changepoints,
    n_nodes = n_nodes, popularity_0_1 = popularity_0_1,
    popularity_1_0 = popularity_1_0, estimate = estimate,
    simultaneous_interactions = simultaneous_interactions,
    verbose = FALSE, it_max = it_max, accelerated = accelerated
  )

  saveRDS(list(res_simulation,
    coef_0_1 = coef_0_1, coef_1_0 = coef_1_0,
    popularity_0_1 = popularity_0_1,
    popularity_1_0 = popularity_1_0,
    baseline_0_1_gt = baseline_0_1_gt,
    baseline_1_0_gt = baseline_1_0_gt
  ), file = paste0(path, "simulation_3_", n_nodes, ".RDS"))
}

simulation_m <- function(x, cutpoints, formula_0_1, formula_1_0, coef_0_1, coef_1_0, popularity_0_1,
                         popularity_1_0, n_nodes, nr = FALSE, simultaneous_interactions, verbose) {
  events <- dem.simulate(
    formula_0_1 = formula_0_1,
    formula_1_0 = formula_1_0,
    coef_0_1 = c(coef_0_1),
    coef_1_0 = c(coef_1_0),
    coef_degree_0_1 = popularity_0_1,
    coef_degree_1_0 = popularity_1_0,
    n_nodes = n_nodes,
    n_events = max(cutpoints),
    verbose = verbose,
    simultaneous_interactions = simultaneous_interactions,
    seed = x
  )


  sec_tmp <- c()
  mem_tmp <- c()
  events_tmp <- c()
  sec_tmp_nr <- c()
  mem_tmp_nr <- c()
  i <- 1
  for (i in 1:length(cutpoints)) {
    cat("Cutpoint ", i, "\n")
    if (nrow(events) >= cutpoints[i]) {
      tmp_events <- events[1:cutpoints[i], ]
      if (i == length(cutpoints)) {
        exogenous_end <- cutpoints[i]
      } else {
        exogenous_end <- events[(cutpoints[i] + 1), 1]
      }
      res_tmp <- peakRAM::peakRAM(est_mm <- dem(
        events = tmp_events,
        formula_0_1 = formula_0_1,
        formula_1_0 = formula_1_0,
        n_nodes = n_nodes,
        exogenous_end = exogenous_end,
        semiparametric = FALSE,
        training_start = 0,
        simultaneous_interactions = simultaneous_interactions,
        control = control.redeem(
          subsample = 1,
          verbose = verbose,
          weighting = FALSE,
          accelerated = FALSE,
          estimate = "Blockwise",
          it_max = 100,
          tol = 1e-20, check_matrix = FALSE
        )
      ))
      if (nr) {
        res_nr <- peakRAM::peakRAM(est_nr <- dem(
          events = tmp_events,
          formula_0_1 = formula_0_1,
          formula_1_0 = formula_1_0,
          n_nodes = n_nodes,
          semiparametric = FALSE,
          training_start = 0,
          exogenous_end = events[(cutpoints[i] + 1), 1],
          simultaneous_interactions = simultaneous_interactions,
          control = control.redeem(
            subsample = 1,
            verbose = verbose,
            weighting = FALSE,
            accelerated = FALSE,
            estimate = "NR",
            it_max = 100,
            tol = 1e-20, check_matrix = FALSE
          )
        ))
        sec_tmp_nr[i] <- res_nr$Elapsed_Time_sec
        mem_tmp_nr[i] <- res_nr$Peak_RAM_Used_MiB
      }
      sec_tmp[i] <- res_tmp$Elapsed_Time_sec
      mem_tmp[i] <- res_tmp$Peak_RAM_Used_MiB
      events_tmp[i] <- nrow(tmp_events)
    }
  }
  if (nr) {
    return(data.frame(
      time_needed = sec_tmp,
      peak_mem = mem_tmp,
      number_events = events_tmp,
      peak_mem_nr = mem_tmp_nr,
      time_needed_nr = sec_tmp_nr
    ))
  } else {
    return(data.frame(
      time_needed = sec_tmp,
      peak_mem = mem_tmp,
      number_events = events_tmp
    ))
  }
}

run_simulation_for_n_node <- function(seed, n_nodes, simultaneous_interactions, it_max, accelerated, estimate, K, verbose, clust, path) {
  set.seed(seed)
  continuous_cov <- rnorm(n = n_nodes)
  continuous_cov <- outer(X = continuous_cov, Y = continuous_cov, FUN = function(x, y) {
    abs(x - y)
  })
  categorical_cov <- sample(x = 1:3, size = n_nodes, replace = T)
  categorical_cov <- outer(X = categorical_cov, Y = categorical_cov, FUN = function(x, y) {
    x == y
  })

  popularity_0_1 <- rnorm(n = n_nodes, sd = 1) - (6 + 0.1 * log(n_nodes))
  popularity_1_0 <- rnorm(n = n_nodes, sd = 1) - (1.6 + 0.1 * log(n_nodes))


  time_changepoints <- seq(0, 10000, length.out = 10)[-c(1, 10)]

  baseline_0_1_gt <- -seq(from = 0, to = -0.1, length.out = length(time_changepoints) + 1)[-1]
  baseline_1_0_gt <- seq(from = 0, to = -0.1, length.out = length(time_changepoints) + 1)[-1]
  formula_0_1 <- ~ degree + current_common_partners(transformation = "log") +
    dyadic_cov(data = continuous_cov) + dyadic_cov(data = categorical_cov) + baseline(changepoints = time_changepoints)
  coef_0_1 <- c(
    "current_common_partners" = -0.5,
    "continuous_cov" = 1,
    "categorical_cov" = 0.5
  )
  formula_1_0 <- ~ degree + inertia(transformation = "log") +
    dyadic_cov(data = continuous_cov) + dyadic_cov(data = categorical_cov) + baseline(changepoints = time_changepoints)
  coef_1_0 <- c(
    "inertia" = 0.5,
    "continuous_cov" = 0.5, "categorical_cov" = 0.5
  )
  # Load the continuous and categorical variable from the working directory to the cluster
  clusterExport(clust,
    varlist = c("continuous_cov", "categorical_cov"),
    envir = environment()
  )
  # debugonce(simulation_run)
  # browser()
  # for(i in 1:100){
  #   paste("Iteration", i,"\n")
  # browser()
  # debugonce(simulation_run)
  # i <- 48
  #   simulation_run(x = i, formula_0_1 = formula_0_1, formula_1_0 = formula_1_0,
  #                  coef_0_1 = coef_0_1, coef_1_0 = coef_1_0,
  #                #                  baseline_0_1_gt = baseline_0_1_gt,
  #                  baseline_1_0_gt = baseline_1_0_gt,
  #                  n_nodes = n_nodes, popularity_0_1 = popularity_0_1,
  #                  popularity_1_0 = popularity_1_0,estimate=estimate,
  #                  simultaneous_interactions = simultaneous_interactions,
  #                  verbose = FALSE, it_max = it_max,accelerated = accelerated)
  # }
  # browser()

  res_simulation <- parLapply(
    cl = clust, X = 1:(20 * K), fun = simulation_run,
    formula_0_1 = formula_0_1, formula_1_0 = formula_1_0,
    coef_0_1 = coef_0_1, coef_1_0 = coef_1_0,
    baseline_0_1_gt = baseline_0_1_gt,
    baseline_1_0_gt = baseline_1_0_gt,
    n_nodes = n_nodes, popularity_0_1 = popularity_0_1,
    popularity_1_0 = popularity_1_0, estimate = estimate,
    simultaneous_interactions = simultaneous_interactions,
    verbose = FALSE, it_max = it_max, accelerated = accelerated
  )

  saveRDS(
    list(res_simulation,
      coef_0_1 = coef_0_1,
      coef_1_0 = coef_1_0,
      popularity_0_1 = popularity_0_1,
      popularity_1_0 = popularity_1_0,
      baseline_0_1_gt = baseline_0_1_gt,
      baseline_1_0_gt = baseline_1_0_gt
    ),
    file = paste0(path, "simulation_1_", n_nodes, ".RDS")
  )
}


get_error <- function(res, coef_0_1,
                      popularity_0_1,
                      baseline_0_1_gt,
                      coef_1_0,
                      popularity_1_0,
                      baseline_1_0_gt) {
  res_core_0_1 <- do.call(
    "rbind",
    lapply(res, function(x) {
      data.frame(
        t(x$coef_0_1_core)
      )
    })
  )
  # Set unidentifable coefs to NA
  res_core_0_1[res_core_0_1 < -20] <- NA

  res_core_1_0 <- do.call(
    "rbind",
    lapply(res, function(x) {
      data.frame(
        t(x$coef_1_0_core)
      )
    })
  )
  res_core_1_0[res_core_1_0 < -20] <- NA
  res_degree_0_1 <- do.call("rbind", lapply(res, function(x) {
    res <- data.frame(t(x$coef_0_1_degree))
    res[x$unidentifiable_0_1] <- NA
    res
  }))
  # plot(popularity_0_1,colMeans(res_degree_0_1))

  res_degree_1_0 <- do.call("rbind", lapply(res, function(x) {
    res <- data.frame(t(x$coef_1_0_degree))
    res[x$unidentifiable_1_0] <- NA
    res
  }))
  # browser()
  res_time_0_1 <- rbindlist(
    lapply(lapply(res, function(x) {
      data.frame(
        t(x$coef_0_1_time)
      )
    }), as.list), # each vector → one-row named list
    fill = TRUE, use.names = TRUE
  )
  res_time_1_0 <- rbindlist(
    lapply(lapply(res, function(x) {
      data.frame(
        t(x$coef_1_0_time)
      )
    }), as.list), # each vector → one-row named list
    fill = TRUE, use.names = TRUE
  )

  coefs_0_1 <- cbind(res_core_0_1, res_degree_0_1, res_time_0_1)
  coefs_1_0 <- cbind(res_core_1_0, res_degree_1_0, res_time_1_0)

  abs_diff_0_1 <- sweep(coefs_0_1, 2, c(
    coef_0_1, popularity_0_1,
    baseline_0_1_gt
  ), FUN = "-")
  abs_diff_1_0 <- sweep(coefs_1_0, 2, c(
    coef_1_0, popularity_1_0,
    baseline_1_0_gt
  ), FUN = "-")

  error_inf_0_1 <- apply(abs_diff_0_1, MARGIN = 1, FUN = function(x) norm(as.matrix(x[!is.na(x)]), type = "I"))
  error_inf_1_0 <- apply(abs_diff_1_0, MARGIN = 1, FUN = function(x) norm(as.matrix(x[!is.na(x)]), type = "I"))
  error_rmse_0_1 <- apply(abs_diff_0_1, MARGIN = 1, FUN = function(x) sqrt(mean(x[!is.na(x)]^2)))
  error_rmse_1_0 <- apply(abs_diff_1_0, MARGIN = 1, FUN = function(x) sqrt(mean(x[!is.na(x)]^2)))

  abs_diff_core_0_1 <- (sweep(res_core_0_1, 2, c(coef_0_1), FUN = "-"))
  error_inf_core_0_1 <- apply(abs_diff_core_0_1, MARGIN = 2, FUN = function(x) max(abs(x[!is.na(x)])))
  error_rmse_core_0_1 <- apply(abs_diff_core_0_1, MARGIN = 2, FUN = function(x) sqrt(mean(x[!is.na(x)]^2)))

  abs_diff_core_1_0 <- (sweep(res_core_1_0, 2, c(coef_1_0), FUN = "-"))

  error_inf_core_1_0 <- apply(abs_diff_core_1_0, MARGIN = 2, FUN = function(x) max(abs(x[!is.na(x)])))
  error_rmse_core_1_0 <- apply(abs_diff_core_1_0, MARGIN = 2, FUN = function(x) sqrt(mean(x[!is.na(x)]^2)))


  error_rmse_alt_core_0_1 <- apply(abs_diff_core_0_1, MARGIN = 1, FUN = function(x) sqrt(mean(x[!is.na(x)]^2, na.rm = TRUE)))
  error_rmse_alt_core_1_0 <- apply(abs_diff_core_1_0, MARGIN = 1, FUN = function(x) sqrt(mean(x[!is.na(x)]^2, na.rm = TRUE)))

  diff_degree_0_1 <- (sweep(as.matrix(res_degree_0_1), 2, c(popularity_0_1), FUN = "-"))
  error_rmse_degree_0_1 <- apply(diff_degree_0_1, MARGIN = 2, FUN = function(x) sqrt(mean(x[!is.na(x)]^2, na.rm = TRUE)))
  error_inf_degree_0_1 <- apply(diff_degree_0_1, MARGIN = 2, FUN = function(x) max(abs(x[!is.na(x)])))
  error_rmse_alt_degree_0_1 <- apply(diff_degree_0_1, MARGIN = 1, FUN = function(x) sqrt(mean(x[!is.na(x)]^2, na.rm = TRUE)))
  error_inf_alt_degree_0_1 <- apply(diff_degree_0_1, MARGIN = 1, FUN = function(x) max(abs(x[!is.na(x)])))

  diff_degree_1_0 <- (sweep(as.matrix(res_degree_1_0), 2, c(popularity_1_0), FUN = "-"))
  error_rmse_degree_1_0 <- apply(diff_degree_1_0, MARGIN = 2, FUN = function(x) sqrt(mean(x[!is.na(x)]^2, na.rm = TRUE)))
  error_inf_degree_1_0 <- apply(diff_degree_1_0, MARGIN = 2, FUN = function(x) max(abs(x[!is.na(x)])))

  error_rmse_alt_degree_1_0 <- apply(diff_degree_1_0, MARGIN = 1, FUN = function(x) sqrt(mean(x[!is.na(x)]^2, na.rm = TRUE)))
  error_inf_alt_degree_1_0 <- apply(diff_degree_1_0, MARGIN = 1, FUN = function(x) max(abs(x[!is.na(x)])))


  diff_time_0_1 <- (sweep(as.matrix(res_time_0_1), 2, c(baseline_0_1_gt), FUN = "-"))
  error_rmse_time_0_1 <- apply(diff_time_0_1, MARGIN = 2, FUN = function(x) sqrt(mean(x[!is.na(x)]^2, na.rm = TRUE)))
  error_inf_time_0_1 <- apply(diff_time_0_1, MARGIN = 2, FUN = function(x) max(abs(x[!is.na(x)])))
  diff_time_1_0 <- (sweep(as.matrix(res_time_1_0), 2, c(baseline_1_0_gt), FUN = "-"))
  error_rmse_time_1_0 <- apply(diff_time_1_0, MARGIN = 2, FUN = function(x) sqrt(mean(x[!is.na(x)]^2, na.rm = TRUE)))
  error_inf_time_1_0 <- apply(diff_time_1_0, MARGIN = 2, FUN = function(x) max(abs(x[!is.na(x)])))
  error_rmse_alt_time_1_0 <- apply(diff_time_1_0, MARGIN = 1, FUN = function(x) sqrt(mean(x[!is.na(x)]^2, na.rm = TRUE)))
  error_rmse_alt_time_0_1 <- apply(diff_time_0_1, MARGIN = 1, FUN = function(x) sqrt(mean(x[!is.na(x)]^2, na.rm = TRUE)))


  # plot(error_rmse_time_1_0,baseline_1_0_gt)
  # mean(error_rmse_degree_1_0)
  #
  # plot(error_rmse_degree_0_1,error_rmse_degree_1_0 )
  # plot(error_rmse_degree_1_0,popularity_1_0)
  # plot(error_rmse_degree_1_0,popularity_0_1 )

  # mean(diff_degree_0_1, na.rm = TRUE)
  # colMeans(abs_diff_core_0_1)
  # colMeans(abs_diff_core_1_0)
  #


  cp_array_0_1 <- do.call(rbind, lapply(res, function(x) data.frame(t(x$cp_0_1))))
  cp_array_1_0 <- do.call(rbind, lapply(res, function(x) data.frame(t(x$cp_1_0))))

  return(list(
    error_inf_0_1 = error_inf_0_1, error_inf_1_0 = error_inf_1_0,
    error_rmse_0_1 = error_rmse_0_1, error_rmse_1_0 = error_rmse_1_0,
    error_inf_core_0_1 = error_inf_core_0_1, error_rmse_core_0_1 = error_rmse_core_0_1,
    error_inf_core_1_0 = error_inf_core_1_0, error_rmse_core_1_0 = error_rmse_core_1_0,
    abs_diff_core_0_1 = abs_diff_core_0_1,
    abs_diff_core_1_0 = abs_diff_core_1_0,
    error_rmse_alt_core_0_1 = error_rmse_alt_core_0_1,
    error_rmse_alt_core_1_0 = error_rmse_alt_core_1_0,
    error_rmse_alt_time_1_0 = error_rmse_alt_time_1_0,
    error_rmse_alt_time_0_1 = error_rmse_alt_time_0_1,
    cp_0_1 = colMeans(cp_array_0_1, na.rm = TRUE),
    cp_1_0 = colMeans(cp_array_1_0, na.rm = TRUE),
    error_rmse_degree_0_1 = error_rmse_degree_0_1,
    error_inf_degree_0_1 = error_inf_degree_0_1,
    error_rmse_alt_degree_0_1 = error_rmse_alt_degree_0_1,
    error_inf_alt_degree_0_1 = error_inf_alt_degree_0_1,
    error_rmse_degree_1_0 = error_rmse_degree_1_0,
    error_inf_degree_1_0 = error_inf_degree_1_0,
    error_rmse_time_0_1 = error_rmse_time_0_1,
    error_inf_time_0_1 = error_inf_time_0_1,
    error_rmse_time_1_0 = error_rmse_time_1_0,
    error_inf_time_1_0 = error_inf_time_1_0,
    error_rmse_alt_degree_1_0 = error_rmse_alt_degree_1_0,
    error_inf_alt_degree_1_0 = error_inf_alt_degree_1_0,
    ave_1_0 = colMeans(res_core_1_0, na.rm = TRUE),
    ave_0_1 = colMeans(res_core_0_1, na.rm = TRUE)
  ))
}


simulation_run_simple <- function(x, formula_0_1, formula_1_0, coef_0_1, coef_1_0,
                                  time_changepoints, baseline_0_1_gt, baseline_1_0_gt,
                                  n_nodes, popularity_0_1, popularity_1_0,
                                  simultaneous_interactions, verbose,
                                  it_max = 1000, accelerated = FALSE,
                                  estimate = "NR") {
  set.seed(x)
  cat(x, " Iteration \n")
  # verbose <- TRUE time = 10000,

  events <- dem.simulate(
    formula_0_1 = formula_0_1,
    formula_1_0 = formula_1_0,
    coef_0_1 = c(coef_0_1),
    coef_1_0 = c(coef_1_0), n_events = 20000,
    coef_degree_0_1 = popularity_0_1,
    coef_degree_1_0 = popularity_1_0,
    n_nodes = n_nodes,
    verbose = verbose,
    simultaneous_interactions = simultaneous_interactions
  )
  # plot(events[,1])
  # plot(events[events[,4]==1,1])
  # plot(events[events[,4]==0,1])
  # debugonce(estimate_mm)
  models <- dem(
    events = events,
    formula_0_1 = formula_0_1,
    formula_1_0 = formula_1_0,
    n_nodes = n_nodes,
    semiparametric = FALSE,
    training_start = 0,
    simultaneous_interactions = simultaneous_interactions,
    control = control.redeem(
      subsample = 1,
      verbose = verbose,
      weighting = FALSE,
      accelerated = accelerated,
      estimate = "Blockwise",
      it_max = it_max,
      tol = 0.0001, check_matrix = FALSE
    )
  )
  unidentifiable_0_1 <- as.vector(which(table(factor(as.vector(events[events[, 4] == 1, c(2, 3)]), 1:(n_nodes))) == 0))
  unidentifiable_1_0 <- as.vector(which(table(factor(as.vector(events[events[, 4] == 0, c(2, 3)]), 1:(n_nodes))) == 0))

  # summary(models)
  # summary(models_alt)
  return(list(
    unidentifiable_0_1 = unidentifiable_0_1,
    unidentifiable_1_0 = unidentifiable_1_0,
    cp_1_0 = (coef_1_0 < models$model_1_0$est_core + qnorm(0.975) * sqrt(diag(models$model_1_0$covariance))) +
      (coef_1_0 > models$model_1_0$est_core - qnorm(0.975) * sqrt(diag(models$model_1_0$covariance))) == 2,
    cp_0_1 = (coef_0_1 < models$model_0_1$est_core + qnorm(0.975) * sqrt(diag(models$model_0_1$covariance))) +
      (coef_0_1 > models$model_0_1$est_core - qnorm(0.975) * sqrt(diag(models$model_0_1$covariance))) == 2,
    coef_1_0_core = models$model_1_0$est_core,
    coef_0_1_core = models$model_0_1$est_core,
    coef_1_0_degree = models$model_1_0$est_degree,
    coef_0_1_degree = models$model_0_1$est_degree,
    coef_1_0_time = models$model_1_0$est_time,
    coef_0_1_time = models$model_0_1$est_time
  ))
  # cp_1_0_newton <- (coef_1_0 < models$model_1_0$coefficient[1:length(coef_1_0)] +
  #                     qnorm(0.975)*sqrt(diag(summary(models$model_1_0)$cov.scaled[1:length(coef_1_0),1:length(coef_1_0)]))) +
  #   (coef_1_0 > models$model_1_0$coefficient[1:length(coef_1_0)] -
  #      qnorm(0.975)*sqrt(diag(summary(models$model_1_0)$cov.scaled[1:length(coef_1_0),1:length(coef_1_0)]))) ==2
  # cp_0_1_newton <- (coef_0_1 < models$model_0_1$coefficient[1:length(coef_0_1)] +
  #                     qnorm(0.975)*sqrt(diag(summary(models$model_0_1)$cov.scaled[1:length(coef_0_1),1:length(coef_0_1)]))) +
  #   (coef_0_1 > models$model_0_1$coefficient[1:length(coef_0_1)] -
  #      qnorm(0.975)*sqrt(diag(summary(models$model_0_1)$cov.scaled[1:length(coef_0_1),1:length(coef_0_1)]))) ==2
  #
  # cat(x," Finished \n")
  #
  # return(list(cp_0_1 = cp_0_1_newton,
  #             cp_1_0 = cp_1_0_newton,
  #             coef_0_1_core = models$model_0_1$coefficients[1:length(coef_0_1)],
  #             coef_1_0_core = models$model_1_0$coefficients[1:length(coef_1_0)]))
}

simulation_model_selection <- function(x, formula_0_1, formula_1_0, coef_0_1, coef_1_0,
                                       time_changepoints, baseline_0_1_gt, baseline_1_0_gt,
                                       n_nodes, popularity_0_1, popularity_1_0,
                                       simultaneous_interactions, verbose,
                                       it_max = 1000, accelerated = FALSE,
                                       estimate = "NR") {
  set.seed(x)
  cat(x, " Iteration \n")
  # verbose <- TRUE
  # debugonce(dem.simulate)
  events <- dem.simulate(
    formula_0_1 = update(formula_0_1, . ~ . + baseline(changepoints = time_changepoints)),
    formula_1_0 = update(formula_1_0, . ~ . + baseline(changepoints = time_changepoints)),
    coef_0_1 = c(coef_0_1),
    coef_1_0 = c(coef_1_0),
    coef_degree_0_1 = popularity_0_1,
    coef_degree_1_0 = popularity_1_0,
    n_nodes = n_nodes,
    time = 10000,
    baseline_0_1 = baseline_0_1_gt,
    baseline_1_0 = baseline_1_0_gt,
    verbose = verbose,
    simultaneous_interactions = simultaneous_interactions,
    seed = x
  )
  # saveRDS(events, file = "../events_tmp.RDS")
  #
  # dim(events)
  # max(events[, 1])

  # cor(models$data_0_1$event, models$data_0_1$prediction)
  # cor(models_tmp$data_0_1$event, models_tmp$data_0_1$prediction)
  # eval_llh_pois(outcome = models$data_0_1$event, mean = models$data_0_1$prediction)
  # eval_llh_pois(outcome = models_tmp$data_0_1$event, mean = models_tmp$data_0_1$prediction)
  # calc_llh_scaled(pred = models$data_0_1$prediction,
  #                       delta = models$data_0_1$event,
  #                       pair_id =  models$data_0_1$pair_id)
  #
  # calc_llh_scaled(pred = models_tmp$data_0_1$prediction,
  #                 delta = models_tmp$data_0_1$event,
  #                 pair_id =  models_tmp$data_0_1$pair_id)
  #
  models <- dem(
    events = events,
    formula_0_1 = update(formula_0_1, . ~ . + baseline(changepoints = time_changepoints)),
    formula_1_0 = update(formula_1_0, . ~ . + baseline(changepoints = time_changepoints)),
    n_nodes = n_nodes,
    exogenous_end = 10000,
    semiparametric = FALSE,
    training_start = 0,
    simultaneous_interactions = simultaneous_interactions,
    control = control.redeem(
      subsample = 1,
      verbose = verbose,
      weighting = FALSE,
      accelerated = FALSE,
      estimate = "Blockwise",
      it_max = it_max,
      tol = 0.00001, check_matrix = FALSE
    )
  )

  # plot(popularity_0_1, models$model_0_1$est_degree)
  terms_0_1 <- terms(formula_0_1)
  terms_1_0 <- terms(formula_1_0)
  n_terms_0_1 <- length(attr(terms_0_1, "term.labels"))
  n_terms_1_0 <- length(attr(terms_1_0, "term.labels"))

  info_model <- summary(models)
  BIC <- numeric(n_terms_1_0 + n_terms_0_1)
  BIC[1] <- info_model$BIC
  AIC <- numeric(n_terms_1_0 + n_terms_0_1)
  AIC[1] <- info_model$AIC

  for (i in 2:(n_terms_1_0 + n_terms_0_1)) {
    cat("Model selection iteration ", i, "\n")
    # if(i == 7){
    #   browser()
    # }
    tmp_drop <- i:(n_terms_1_0 + n_terms_0_1)
    tmp_drop_0_1 <- tmp_drop[tmp_drop <= n_terms_0_1]
    tmp_drop_1_0 <- tmp_drop[tmp_drop > n_terms_0_1] - n_terms_0_1
    formula_0_1_tmp <- formula(drop.terms(terms_0_1, dropx = tmp_drop_0_1))
    if (length(tmp_drop_1_0) == n_terms_1_0) {
      formula_1_0_tmp <- ~degree
    } else {
      formula_1_0_tmp <- formula(drop.terms(terms_1_0, dropx = tmp_drop_1_0))
    }
    models_tmp <- dem(
      events = events,
      formula_0_1 = update(formula_0_1_tmp, . ~ . + baseline(changepoints = time_changepoints)),
      formula_1_0 = update(formula_1_0_tmp, . ~ . + baseline(changepoints = time_changepoints)),
      n_nodes = n_nodes,
      exogenous_end = 10000,
      semiparametric = FALSE,
      training_start = 0,
      simultaneous_interactions = simultaneous_interactions,
      control = control.redeem(
        subsample = 1,
        verbose = verbose,
        weighting = FALSE,
        accelerated = FALSE,
        estimate = "Blockwise",
        it_max = it_max,
        tol = 0.001, check_matrix = FALSE
      )
    )
    info_model_tmp <- summary(models_tmp)

    BIC[i] <- info_model_tmp$BIC
    AIC[i] <- info_model_tmp$AIC
  }


  cat(x, " Finished \n")

  # summary(models)
  unidentifiable_0_1 <- as.vector(which(table(factor(as.vector(events[events[, 4] == 1, c(2, 3)]), 1:(n_nodes))) == 0))
  unidentifiable_1_0 <- as.vector(which(table(factor(as.vector(events[events[, 4] == 0, c(2, 3)]), 1:(n_nodes))) == 0))
  cat(x, " Finished \n")

  return(list(
    unidentifiable_0_1 = unidentifiable_0_1, unidentifiable_1_0 = unidentifiable_1_0,
    cp_1_0 = (coef_1_0 < models$model_1_0$est_core + qnorm(0.975) * sqrt(diag(models$model_1_0$covariance))) +
      (coef_1_0 > models$model_1_0$est_core - qnorm(0.975) * sqrt(diag(models$model_1_0$covariance))) == 2,
    cp_0_1 = (coef_0_1 < models$model_0_1$est_core + qnorm(0.975) * sqrt(diag(models$model_0_1$covariance))) +
      (coef_0_1 > models$model_0_1$est_core - qnorm(0.975) * sqrt(diag(models$model_0_1$covariance))) == 2,
    covariance_0_1 = models$model_0_1$covariance,
    covariance_1_0 = models$model_1_0$covariance,
    coef_1_0_core = models$model_1_0$est_core,
    coef_0_1_core = models$model_0_1$est_core,
    coef_1_0_degree = models$model_1_0$est_degree,
    coef_0_1_degree = models$model_0_1$est_degree,
    coef_1_0_time = models$model_1_0$est_time,
    coef_0_1_time = models$model_0_1$est_time,
    BIC = BIC, AIC = AIC,
    coef_0_1 = coef_0_1, coef_1_0 = coef_1_0,
    baseline_0_1_gt = baseline_0_1_gt, baseline_1_0_gt = baseline_1_0_gt,
    popularity_0_1 = popularity_0_1, popularity_1_0 = popularity_1_0
  ))
}


simulation_baseline <- function(x, formula_0_1, formula_1_0, coef_0_1, coef_1_0,
                                time_changepoints,
                                baseline_0_1_gt, baseline_1_0_gt,
                                n_nodes, popularity_0_1, popularity_1_0,
                                simultaneous_interactions, verbose,
                                it_max = 1000,
                                accelerated = FALSE, grid_cells) {
  set.seed(x)
  cat(x, " Iteration \n")
  now <- Sys.time()
  events <- dem.simulate(
    formula_0_1 = update(formula_0_1, . ~ . + baseline(changepoints = time_changepoints)),
    formula_1_0 = update(formula_1_0, . ~ . + baseline(changepoints = time_changepoints)),
    coef_0_1 = c(coef_0_1),
    coef_1_0 = c(coef_1_0),
    coef_degree_0_1 = popularity_0_1,
    coef_degree_1_0 = popularity_1_0,
    n_nodes = n_nodes, time = 10000,
    baseline_0_1 = baseline_0_1_gt,
    baseline_1_0 = baseline_1_0_gt,
    verbose = F,
    simultaneous_interactions = simultaneous_interactions,
    seed = x
  )
  passed <- Sys.time() - now
  cat("Simulation done in ", passed, "\n")

  unidentifiable_0_1 <- as.vector(which(table(factor(as.vector(events[events[, 4] == 1, c(2, 3)]), 1:(n_nodes))) == 0))
  unidentifiable_1_0 <- as.vector(which(table(factor(as.vector(events[events[, 4] == 0, c(2, 3)]), 1:(n_nodes))) == 0))
  # i <- 500
  cp_0_1_list <- list()
  cp_1_0_list <- list()
  est_0_1_list <- list()
  est_1_0_list <- list()
  needed_time <- list()
  l <- 1
  for (i in grid_cells) {
    cat(l, "\n")
    time_changepoints_tmp <- seq(0, 10000, length.out = i + 1)[-1]
    now <- Sys.time()
    cp_vals <- (time_changepoints_tmp)[-length(time_changepoints_tmp)]
    formula_0_1_tmp <- update(formula_0_1, . ~ . + baseline(changepoints = cp_vals))
    environment(formula_0_1_tmp) <- environment()
    formula_1_0_tmp <- update(formula_1_0, . ~ . + baseline(changepoints = cp_vals))
    environment(formula_1_0_tmp) <- environment()
    now <- Sys.time()
    models_tmp <- dem(
      events = events,
      formula_0_1 = formula_0_1_tmp,
      formula_1_0 = formula_1_0_tmp,
      n_nodes = n_nodes,
      exogenous_end = 10000,
      semiparametric = FALSE,
      training_start = 0,
      simultaneous_interactions = simultaneous_interactions,
      control = control.redeem(
        subsample = 1,
        verbose = verbose,
        weighting = FALSE,
        accelerated = FALSE,
        estimate = "Blockwise",
        it_max = c(it_max, 2*it_max) + i,
        tol = 0.001, check_matrix = FALSE
      )
    )

    passed <- Sys.time() - now
    cat("Estimation done in ", passed, "\n")

    needed_time[[l]] <- Sys.time() - now
    units(needed_time[[l]]) <- "mins"
    needed_time[[l]] <- needed_time[[l]] / (nrow(models_tmp$model_0_1$est_degree_hist) + nrow(models_tmp$model_1_0$est_degree_hist))
    est_0_1_list[[l]] <- models_tmp$model_0_1$est_core
    est_1_0_list[[l]] <- models_tmp$model_1_0$est_core
    cp_1_0_list[[l]] <- (coef_1_0 < models_tmp$model_1_0$est_core + qnorm(0.975) * sqrt(diag(models_tmp$model_1_0$covariance))) +
      (coef_1_0 > models_tmp$model_1_0$est_core - qnorm(0.975) * sqrt(diag(models_tmp$model_1_0$covariance))) == 2
    cp_0_1_list[[l]] <- (coef_0_1 < models_tmp$model_0_1$est_core + qnorm(0.975) * sqrt(diag(models_tmp$model_0_1$covariance))) +
      (coef_0_1 > models_tmp$model_0_1$est_core - qnorm(0.975) * sqrt(diag(models_tmp$model_0_1$covariance))) == 2
    l <- l + 1
  }


  cat(x, " Finished \n")
  return(list(
    n_events_0_1 = sum(events[, 4] == 1),
    n_events_1_0 = sum(events[, 4] == 0),
    cp_0_1_list = cp_0_1_list,
    cp_1_0_list = cp_1_0_list,
    needed_time = needed_time,
    est_0_1_list = est_0_1_list,
    est_1_0_list = est_1_0_list,
    grid_cells = grid_cells,
    unidentifiable_0_1 = unidentifiable_0_1,
    unidentifiable_1_0 = unidentifiable_1_0,
    coef_0_1 = coef_0_1, coef_1_0 = coef_1_0,
    baseline_0_1_gt = baseline_0_1_gt,
    baseline_1_0_gt = baseline_1_0_gt,
    popularity_0_1 = popularity_0_1,
    popularity_1_0 = popularity_1_0
  ))
}


simulation_run <- function(x, formula_0_1, formula_1_0, coef_0_1, coef_1_0,
                           time_changepoints, baseline_0_1_gt, baseline_1_0_gt,
                           n_nodes, popularity_0_1, popularity_1_0,
                           simultaneous_interactions, verbose,
                           it_max = 1000, accelerated = FALSE,
                           estimate = TRUE) {
  set.seed(x)
  cat(x, " Iteration \n")
  # verbose <- TRUE
  # debugonce(dem.simulate)
  events <- dem.simulate(
    formula_0_1 = update(formula_0_1, . ~ . + baseline(changepoints = time_changepoints)),
    formula_1_0 = update(formula_1_0, . ~ . + baseline(changepoints = time_changepoints)),
    coef_0_1 = c(coef_0_1),
    coef_1_0 = c(coef_1_0),
    coef_degree_0_1 = popularity_0_1,
    coef_degree_1_0 = popularity_1_0,
    n_nodes = n_nodes, time = 10000,
    baseline_0_1 = baseline_0_1_gt,
    baseline_1_0 = baseline_1_0_gt,
    verbose = verbose,
    simultaneous_interactions = simultaneous_interactions,
    seed = x
  )

  # summary(models)
  # plot(events[,1])
  if (estimate == "NR") {
    model_newton <- dem(
      events = events,
      formula_0_1 = update(formula_0_1, . ~ . + baseline(changepoints = time_changepoints[-length(time_changepoints)], labels = paste0("time_", time_changepoints[-length(time_changepoints)]))),
      formula_1_0 = update(formula_1_0, . ~ . + baseline(changepoints = time_changepoints[-length(time_changepoints)], labels = paste0("time_", time_changepoints[-length(time_changepoints)]))),
      n_nodes = n_nodes,
      exogenous_end = 10000,
      semiparametric = FALSE,
      training_start = 0,
      simultaneous_interactions = simultaneous_interactions,
      control = control.redeem(
        subsample = 1,
        verbose = verbose,
        weighting = FALSE,
        accelerated = accelerated,
        estimate = "NR",
        it_max = it_max,
        tol = 0.0001, check_matrix = FALSE
      )
    )

    cp_1_0_newton <- (coef_1_0 < model_newton$model_1_0$coefficient[1:length(coef_1_0)] +
      qnorm(0.975) * sqrt(diag(summary(model_newton$model_1_0)$cov.scaled[1:length(coef_1_0), 1:length(coef_1_0)]))) +
      (coef_1_0 > model_newton$model_1_0$coefficient[1:length(coef_1_0)] -
        qnorm(0.975) * sqrt(diag(summary(model_newton$model_1_0)$cov.scaled[1:length(coef_1_0), 1:length(coef_1_0)]))) == 2
    cp_0_1_newton <- (coef_0_1 < model_newton$model_0_1$coefficient[1:length(coef_0_1)] +
      qnorm(0.975) * sqrt(diag(summary(model_newton$model_0_1)$cov.scaled[1:length(coef_0_1), 1:length(coef_0_1)]))) +
      (coef_0_1 > model_newton$model_0_1$coefficient[1:length(coef_0_1)] -
        qnorm(0.975) * sqrt(diag(summary(model_newton$model_0_1)$cov.scaled[1:length(coef_0_1), 1:length(coef_0_1)]))) == 2


    newton_info <- list(
      cp_0_1 = as.vector(cp_0_1_newton),
      cp_1_0 = as.vector(cp_1_0_newton),
      coef_0_1_core = as.vector(model_newton$model_0_1$coefficients[1:length(coef_0_1)]),
      coef_1_0_core = as.vector(model_newton$model_1_0$coefficients[1:length(coef_1_0)])
    )
    rm(model_newton)
    gc()
  } else {
    newton_info <- NA
  }

  models <- dem(
    events = events,
    formula_0_1 = update(formula_0_1, . ~ . + baseline(changepoints = (time_changepoints)[-length(time_changepoints)], labels = paste0("time_", (time_changepoints)[-length(time_changepoints)]))),
    formula_1_0 = update(formula_1_0, . ~ . + baseline(changepoints = (time_changepoints)[-length(time_changepoints)], labels = paste0("time_", (time_changepoints)[-length(time_changepoints)]))),
    n_nodes = n_nodes,
    exogenous_end = 10000,
    semiparametric = FALSE,
    training_start = 0,
    simultaneous_interactions = simultaneous_interactions,
    control = control.redeem(
      subsample = 1,
      verbose = verbose,
      weighting = FALSE,
      accelerated = FALSE,
      estimate = "Blockwise",
      it_max = it_max,
      tol = 0.001, check_matrix = FALSE
    )
  )

  # summary(models)
  unidentifiable_0_1 <- as.vector(which(table(factor(as.vector(events[events[, 4] == 1, c(2, 3)]), 1:(n_nodes))) == 0))
  unidentifiable_1_0 <- as.vector(which(table(factor(as.vector(events[events[, 4] == 0, c(2, 3)]), 1:(n_nodes))) == 0))
  cat(x, " Finished \n")

  return(list(
    unidentifiable_0_1 = as.vector(unidentifiable_0_1),
    unidentifiable_1_0 = as.vector(unidentifiable_1_0),
    cp_1_0 = as.vector((coef_1_0 < models$model_1_0$est_core + qnorm(0.975) * sqrt(diag(models$model_1_0$covariance))) +
      (coef_1_0 > models$model_1_0$est_core - qnorm(0.975) * sqrt(diag(models$model_1_0$covariance))) == 2),
    cp_0_1 = as.vector((coef_0_1 < models$model_0_1$est_core + qnorm(0.975) * sqrt(diag(models$model_0_1$covariance))) +
      (coef_0_1 > models$model_0_1$est_core - qnorm(0.975) * sqrt(diag(models$model_0_1$covariance))) == 2),
    coef_1_0_core = as.vector(models$model_1_0$est_core),
    coef_0_1_core = as.vector(models$model_0_1$est_core),
    coef_1_0_degree = as.vector(models$model_1_0$est_degree),
    coef_0_1_degree = as.vector(models$model_0_1$est_degree),
    coef_1_0_time = as.vector(models$model_1_0$est_time),
    coef_0_1_time = as.vector(models$model_0_1$est_time),
    newton_info = newton_info, number_events = nrow(events)
  ))
}


simulation_compare_time_memory <- function(x, formula_0_1, formula_1_0, coef_0_1, coef_1_0,
                                           time_changepoints, baseline_0_1_gt, baseline_1_0_gt,
                                           n_nodes, popularity_0_1, popularity_1_0,
                                           simultaneous_interactions, verbose,
                                           it_max = 1000, accelerated = FALSE,
                                           estimate = "NR") {
  set.seed(x)
  cat(x, " Iteration \n")
  # verbose <- TRUE
  # debugonce(dem.simulate)
  events <- dem.simulate(
    formula_0_1 = update(formula_0_1, . ~ . + baseline(changepoints = time_changepoints)),
    formula_1_0 = update(formula_1_0, . ~ . + baseline(changepoints = time_changepoints)),
    coef_0_1 = c(coef_0_1),
    coef_1_0 = c(coef_1_0),
    coef_degree_0_1 = popularity_0_1,
    coef_degree_1_0 = popularity_1_0,
    n_nodes = n_nodes, time = 10000,
    baseline_0_1 = baseline_0_1_gt,
    baseline_1_0 = baseline_1_0_gt,
    verbose = verbose,
    simultaneous_interactions = simultaneous_interactions,
    seed = x
  )

  # debugonce(estimate_mmt)
  res_mm <- peakRAM::peakRAM(est_mm <- dem(
    events = events,
    formula_0_1 = update(formula_0_1, . ~ . + baseline(changepoints = time_changepoints)),
    formula_1_0 = update(formula_1_0, . ~ . + baseline(changepoints = time_changepoints)),
    n_nodes = n_nodes,
    exogenous_end = 10000,
    semiparametric = FALSE,
    training_start = 0,
    simultaneous_interactions = simultaneous_interactions,
    control = control.redeem(
      subsample = 1,
      verbose = verbose,
      weighting = FALSE,
      accelerated = accelerated,
      estimate = "Blockwise",
      it_max = it_max,
      tol = 0.001, check_matrix = FALSE
    )
  ))

  res_nr <- peakRAM::peakRAM(est_nr <- dem(
    events = events,
    formula_0_1 = update(formula_0_1, . ~ . + baseline(changepoints = time_changepoints)),
    formula_1_0 = update(formula_1_0, . ~ . + baseline(changepoints = time_changepoints)),
    n_nodes = n_nodes,
    exogenous_end = 10000,
    semiparametric = FALSE,
    training_start = 0,
    simultaneous_interactions = simultaneous_interactions,
    control = control.redeem(
      subsample = 1,
      verbose = verbose,
      weighting = FALSE,
      accelerated = accelerated,
      estimate = "NR",
      legacy = TRUE,
      it_max = it_max,
      tol = 0.001, check_matrix = FALSE
    )
  ))
  # res_gd <- peakRAM::peakRAM(est_gd <- dem(
  #   events = events,
  #   formula_0_1 = update(formula_0_1, . ~ . + baseline(changepoints =time_changepoints)),
  #   formula_1_0 = update(formula_1_0, . ~ . + baseline(changepoints =time_changepoints)),
  #   n_nodes = n_nodes,
  #   exogenous_end = (time_changepoints)[length(time_changepoints)],
  #   semiparametric = FALSE,
  #   training_start = 0,
  #   control = control.redeem(
  #     subsample = 1,
  #     verbose = verbose,
  #     weighting = FALSE,
  #     accelerated = accelerated,
  #     simultaneous_interactions = simultaneous_interactions,
  #     estimate = "GD",
  #     it_max = it_max,
  #     tol = 0.001, check_matrix = FALSE
  #   )
  # ))

  time_mm <- res_mm$Elapsed_Time_sec
  time_nr <- res_nr$Elapsed_Time_sec
  # time_gd <- res_gd$Elapsed_Time_sec

  total_ram_mm <- res_mm$Total_RAM_Used_MiB
  total_ram_nr <- res_nr$Total_RAM_Used_MiB
  # total_ram_gd <- res_gd$Total_RAM_Used_MiB

  peak_ram_mm <- res_mm$Peak_RAM_Used_MiB
  peak_ram_nr <- res_nr$Peak_RAM_Used_MiB
  # peak_ram_gd <- res_gd$Peak_RAM_Used_MiB

  unidentifiable_0_1 <- as.vector(which(table(factor(as.vector(events[events[, 4] == 1, c(2, 3)]), 1:(n_nodes))) == 0))
  unidentifiable_1_0 <- as.vector(which(table(factor(as.vector(events[events[, 4] == 0, c(2, 3)]), 1:(n_nodes))) == 0))


  est_mm_0_1 <- c(est_mm$model_0_1$est_core, est_mm$model_0_1$est_degree, est_mm$model_0_1$est_time)
  est_nr_0_1 <- est_nr$model_0_1$coefficients
  # est_gd_0_1 <- est_gd$model_0_1$beta

  est_mm_1_0 <- c(est_mm$model_1_0$est_core, est_mm$model_1_0$est_degree, est_mm$model_1_0$est_time)
  est_nr_1_0 <- est_nr$model_1_0$coefficients
  # est_gd_1_0 <- est_gd$model_1_0$beta

  llh_mm <- est_mm$model_0_1$llh + est_mm$model_1_0$llh
  llh_nr <- as.vector(logLik(est_nr$model_0_1) + logLik(est_nr$model_1_0))
  # llh_gd <- max(est_gd$model_0_1$llh) + max(est_gd$model_1_0$llh)

  cat(x, " Finished \n")

  return(list(
    time_mm = time_mm,
    time_nr = time_nr,
    # time_gd = time_gd,
    total_ram_mm = total_ram_mm,
    total_ram_nr = total_ram_nr,
    # total_ram_gd = total_ram_gd,
    peak_ram_mm = peak_ram_mm,
    peak_ram_nr = peak_ram_nr,
    # peak_ram_gd = peak_ram_gd,
    est_mm_0_1 = est_mm_0_1,
    est_mm_1_0 = est_mm_1_0,
    est_nr_0_1 = est_nr_0_1,
    est_nr_1_0 = est_nr_1_0,
    # est_gd_0_1 = est_gd_0_1,
    # est_gd_1_0 = est_gd_1_0,
    unidentifiable_0_1 = unidentifiable_0_1,
    unidentifiable_1_0 = unidentifiable_1_0, llh_mm, llh_nr
  ))
}

simulation_compare_time_memory_bench <- function(x, formula_0_1, formula_1_0, coef_0_1, coef_1_0,
                                                 time_changepoints, baseline_0_1_gt, baseline_1_0_gt,
                                                 n_nodes, popularity_0_1, popularity_1_0,
                                                 simultaneous_interactions, verbose,
                                                 it_max = 1000, accelerated = FALSE,
                                                 estimate = "NR") {
  set.seed(x)
  cat(x, " Iteration \n")
  events <- dem.simulate(
    formula_0_1 = update(formula_0_1, . ~ . + baseline(changepoints = time_changepoints)),
    formula_1_0 = update(formula_1_0, . ~ . + baseline(changepoints = time_changepoints)),
    coef_0_1 = c(coef_0_1),
    coef_1_0 = c(coef_1_0),
    coef_degree_0_1 = popularity_0_1,
    coef_degree_1_0 = popularity_1_0,
    n_nodes = n_nodes, time = 10000,
    baseline_0_1 = baseline_0_1_gt,
    baseline_1_0 = baseline_1_0_gt,
    verbose = verbose,
    simultaneous_interactions = simultaneous_interactions,
    seed = x
  )

  res <- bench::mark(
    mm = dem(
      events = events,
      formula_0_1 = update(formula_0_1, . ~ . + baseline(changepoints = (time_changepoints)[-length(time_changepoints)], labels = paste0("time_", (time_changepoints)[-length(time_changepoints)]))),
      formula_1_0 = update(formula_1_0, . ~ . + baseline(changepoints = (time_changepoints)[-length(time_changepoints)], labels = paste0("time_", (time_changepoints)[-length(time_changepoints)]))),
      n_nodes = n_nodes,
      exogenous_end = 10000,
      semiparametric = FALSE,
      training_start = 0,
      simultaneous_interactions = simultaneous_interactions,
      control = control.redeem(
        subsample = 1,
        verbose = verbose,
        weighting = FALSE,
        accelerated = accelerated,
        estimate = "Blockwise",
        it_max = it_max, save_hist = FALSE,
        tol = 0.001, check_matrix = FALSE
      )
    ),
    nr = dem(
      events = events,
      formula_0_1 = update(formula_0_1, . ~ . + baseline(changepoints = time_changepoints)),
      formula_1_0 = update(formula_1_0, . ~ . + baseline(changepoints = time_changepoints)),
      n_nodes = n_nodes,
      exogenous_end = 10000,
      semiparametric = FALSE,
      training_start = 0,
      simultaneous_interactions = simultaneous_interactions,
      control = control.redeem(
        subsample = 1,
        verbose = verbose,
        weighting = FALSE,
        accelerated = accelerated,
        estimate = "NR",
        legacy = TRUE,
        it_max = it_max,
        tol = 0.001, check_matrix = FALSE
      )
    ),
    iterations = 1, check = FALSE
  )

  unidentifiable_0_1 <- as.vector(which(table(factor(as.vector(events[events[, 4] == 1, c(2, 3)]), 1:(n_nodes))) == 0))
  unidentifiable_1_0 <- as.vector(which(table(factor(as.vector(events[events[, 4] == 0, c(2, 3)]), 1:(n_nodes))) == 0))

  cat(x, " Finished \n")

  return(list(
    bench_results = res,
    unidentifiable_0_1 = unidentifiable_0_1,
    unidentifiable_1_0 = unidentifiable_1_0
  ))
}

run_simulation_for_n_node_time_bench <- function(seed, n_nodes, simultaneous_interactions, it_max, accelerated, estimate, K, verbose, clust, path) {
  set.seed(seed)
  continuous_cov <- rnorm(n = n_nodes)
  continuous_cov <- outer(X = continuous_cov, Y = continuous_cov, FUN = function(x, y) {
    abs(x - y)
  })
  categorical_cov <- sample(x = 1:3, size = n_nodes, replace = T)
  categorical_cov <- outer(X = categorical_cov, Y = categorical_cov, FUN = function(x, y) {
    x == y
  })

  popularity_0_1 <- rnorm(n = n_nodes, sd = 1) - (6 + 0.1 * log(n_nodes))
  popularity_1_0 <- rnorm(n = n_nodes, sd = 1) - (1.6 + 0.1 * log(n_nodes))


  time_changepoints <- seq(0, 10000, length.out = 10)[-c(1, 10)]

  baseline_0_1_gt <- -seq(from = 0, to = -0.1, length.out = length(time_changepoints) + 1)[-1]
  baseline_1_0_gt <- seq(from = 0, to = -0.1, length.out = length(time_changepoints) + 1)[-1]
  formula_0_1 <- ~ degree + current_common_partners(transformation = "log") +
    dyadic_cov(data = continuous_cov) + dyadic_cov(data = categorical_cov)
  coef_0_1 <- c(
    "current_common_partners" = -0.5,
    "continuous_cov" = 1,
    "categorical_cov" = 0.5
  )
  formula_1_0 <- ~ degree + inertia(transformation = "log") +
    dyadic_cov(data = continuous_cov) + dyadic_cov(data = categorical_cov)
  coef_1_0 <- c(
    "inertia" = 0.5,
    "continuous_cov" = 0.5, "categorical_cov" = 0.5
  )
  # Load the continuous and categorical variable from the working directory to the cluster
  clusterExport(clust,
    varlist = c("continuous_cov", "categorical_cov"),
    envir = environment()
  )

  res_simulation <- parLapply(
    cl = clust, X = 1:(K), fun = simulation_compare_time_memory_bench,
    formula_0_1 = formula_0_1, formula_1_0 = formula_1_0,
    coef_0_1 = coef_0_1, coef_1_0 = coef_1_0,
    baseline_0_1_gt = baseline_0_1_gt,
    baseline_1_0_gt = baseline_1_0_gt,
    time_changepoints = time_changepoints,
    n_nodes = n_nodes, popularity_0_1 = popularity_0_1,
    popularity_1_0 = popularity_1_0, estimate = estimate,
    simultaneous_interactions = simultaneous_interactions,
    verbose = FALSE, it_max = it_max, accelerated = accelerated
  )

  saveRDS(list(res_simulation,
    coef_0_1 = coef_0_1, coef_1_0 = coef_1_0,
    popularity_0_1 = popularity_0_1,
    popularity_1_0 = popularity_1_0,
    baseline_0_1_gt = baseline_0_1_gt,
    baseline_1_0_gt = baseline_1_0_gt
  ), file = paste0(path, "simulation_3_bench_", n_nodes, ".RDS"))
}
