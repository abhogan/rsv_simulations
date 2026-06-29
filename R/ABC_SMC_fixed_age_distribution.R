RSVsim_ABC_SMC_fixed_age_distribution <- function(target,
                                                  epsilon_matrix,
                                                  summary_fun,
                                                  dist_fun,
                                                  prior_fun,
                                                  n_prior_attempts,
                                                  used_seed_matrix,
                                                  prior_dens_fun,
                                                  particle_low,
                                                  particle_up,
                                                  nparticles,
                                                  ncores = 1,
                                                  fitted_parameter_names,
                                                  fixed_parameter_list,
                                                  times = seq(0, 365.25 * 5, 0.25),
                                                  cohort_step_size = 1 / 12 * 365.25,
                                                  warm_up = 365.25 * 4) {
  nparams <- length(fitted_parameter_names)
  ntargets <- length(target)
  if (!all(is.numeric(target))) {
    stop("RSVsim_ABC_SMC_fixed_age_distribution: target must be numeric")
  }
  if (!all(is.numeric(epsilon_matrix))) {
    stop("RSVsim_ABC_SMC_fixed_age_distribution: epsilon_matrix must be numeric")
  }
  if (any(epsilon_matrix <= 0)) {
    stop("RSVsim_ABC_SMC_fixed_age_distribution: increase all epsilon_matrix values above zero")
  }
  if (ncol(epsilon_matrix) != ntargets) {
    stop("RSVsim_ABC_SMC_fixed_age_distribution: incorrect number of columns in epsilon_matrix")
  }
  if (length(particle_up) != length(particle_low) || length(particle_low) != nparams) {
    stop("RSVsim_ABC_SMC_fixed_age_distribution: incorrect number of lower or upper particle boundaries")
  }
  if (ncol(prior_fun(10)) != nparams || nrow(prior_fun(10)) != 10) {
    stop("RSVsim_ABC_SMC_fixed_age_distribution: prior_fun should return a matrix with the number of columns equal to the number of fitted parameters and the number of rows equal to n_prior_attempts")
  }
  if (length(prior_dens_fun(rep(1, nparams))) != nparams) {
    stop("RSVsim_ABC_SMC_fixed_age_distribution: prior_dens_fun should return a vector of the same length as the number of fitted parameters")
  }
  
  G <- nrow(epsilon_matrix)
  n <- 1
  while_fun_SMC <- function(particle, g, w_old, res_old, nparticles, sigma,
                            n_prior_attempts, n, target, epsilon_matrix,
                            fitted_parameter_names, fixed_parameter_list,
                            times, cohort_step_size, warm_up, prior_fun,
                            prior_dens_fun, particle_low, particle_up,
                            used_seed_matrix, ntargets, nparams) {
    used_seed <- used_seed_matrix[g, particle]
    set.seed(used_seed)
    if (g == 1) {
      fitted_parameters <- prior_fun(n_prior_attempts)
    } else {
      p <- sample(seq(1, nparticles), n_prior_attempts, prob = w_old, replace = TRUE)
      fitted_parameters <- as.data.frame(t(sapply(1:n_prior_attempts, function(a) {
        tmvtnorm::rtmvnorm(1,
                          mean = unlist(as.vector(unname(res_old[p[a], , drop = FALSE]))),
                          sigma = sigma,
                          lower = particle_low,
                          upper = particle_up)
      }, simplify = TRUE)))
      if (nparams == 1) {
        fitted_parameters <- t(unname(fitted_parameters))
      }
    }
    i <- 1
    k <- 1
    while (i <= 1) {
      if (k > n_prior_attempts) {
        stop("RSVsim_ABC_SMC_fixed_age_distribution: no particle accepted - increase n_prior_attempts or the tolerance (epsilon_matrix)")
      }
      parameters <- as.numeric(fitted_parameters[k, ])
      parameters_ODE <- RSVsim::RSVsim_update_parameters(fixed_parameter_list, fitted_parameter_names, parameters)
      p_non_zero <- as.numeric(prod(prior_dens_fun(as.vector(parameters))) > 0)
      if (p_non_zero) {
        m <- 0
        distance <- matrix(nrow = n, ncol = ntargets)
        for (j in 1:n) {
          out <- RSVsim_run_model_fixed_age_distribution(parameters = parameters_ODE,
                                                         times = times,
                                                         cohort_step_size = cohort_step_size,
                                                         warm_up = warm_up)
          target_star <- summary_fun(out)
          rm(out)
          distance[j, ] <- dist_fun(target, target_star)
          if (all(distance[j, ] <= epsilon_matrix[g, ])) {
            m <- m + 1
          }
        }
        if (m > 0) {
          res_new_out <- parameters
          w1 <- prod(prior_dens_fun(as.vector(parameters)))
          if (g == 1) {
            w2 <- 1
          } else {
            w2 <- sum(sapply(1:nparticles, function(a) {
              w_old[a] * tmvtnorm::dtmvnorm(as.vector(unname(unlist(res_new_out))),
                                           mean = as.vector(unname(unlist(res_old[a, ]))),
                                           sigma = sigma,
                                           lower = particle_low,
                                           upper = particle_up)
            }))
          }
          w_new_out <- (m / n) * w1 / w2
          i <- i + 1
        }
      }
      k <- k + 1
    }
    return(list(res_new = res_new_out, w_new = w_new_out, seed = used_seed))
  }
  
  if (ncores > 1) {
    cl <- parallel::makePSOCKcluster(ncores)
    base::on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl,
                            varlist = c("while_fun_SMC",
                                        "RSVsim_run_model_fixed_age_distribution",
                                        "fixed_parameter_list",
                                        "fitted_parameter_names",
                                        "times",
                                        "cohort_step_size",
                                        "warm_up",
                                        "target",
                                        "epsilon_matrix",
                                        "prior_fun",
                                        "prior_dens_fun",
                                        "n_prior_attempts",
                                        "summary_fun",
                                        "dist_fun",
                                        "nparams",
                                        "used_seed_matrix",
                                        "particle_low",
                                        "particle_up",
                                        "ntargets"),
                            envir = environment())
    parallel::clusterEvalQ(cl, {
      library(RSVsim)
      library(tidyr)
      library(dplyr)
      library(tmvtnorm)
      base::invisible(data.table::setDTthreads(1))
    })
  } else {
    cl <- NULL
  }
  
  w_list <- vector(mode = "list", length = G)
  res_list <- vector(mode = "list", length = G)
  sigma_list <- vector(mode = "list", length = G)
  res_old <- matrix(nrow = nparticles, ncol = nparams)
  res_new <- matrix(nrow = nparticles, ncol = nparams)
  w_old <- matrix(nrow = nparticles, ncol = 1)
  w_new <- matrix(nrow = nparticles, ncol = 1)
  sigma <- matrix(nrow = nparams, ncol = nparams)
  
  for (g in 1:G) {
    print(paste("Generation", g, "of", G))
    if (ncores > 1) {
      parallel::clusterExport(cl, varlist = c("w_old", "res_old", "sigma", "g"), envir = environment())
    }
    pbapply::pboptions(type = "timer", char = "-")
    res <- pbapply::pblapply(cl = cl,
                             X = 1:nparticles,
                             FUN = while_fun_SMC,
                             g = g,
                             w_old = w_old,
                             res_old = res_old,
                             nparticles = nparticles,
                             sigma = sigma,
                             n_prior_attempts = n_prior_attempts,
                             n = n,
                             target = target,
                             epsilon_matrix = epsilon_matrix,
                             fitted_parameter_names = fitted_parameter_names,
                             fixed_parameter_list = fixed_parameter_list,
                             times = times,
                             cohort_step_size = cohort_step_size,
                             warm_up = warm_up,
                             prior_fun = prior_fun,
                             prior_dens_fun = prior_dens_fun,
                             particle_low = particle_low,
                             particle_up = particle_up,
                             used_seed_matrix = used_seed_matrix,
                             ntargets = ntargets,
                             nparams = nparams)
    w_new <- sapply(res, "[[", "w_new")
    res_new <- do.call(rbind, lapply(res, "[[", "res_new"))
    sigma <- stats::cov(res_new)
    res_old <- res_new
    w_old <- w_new / sum(w_new)
    w_list[[g]] <- w_old
    res_list[[g]] <- res_old
    sigma_list[[g]] <- sigma
  }
  
  return(list(fitted_parameters = res_list, weights = w_list, sigma = sigma_list))
}
