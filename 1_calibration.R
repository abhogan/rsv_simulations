source("R/set_up_variables.r")

# scenario name
sn <- "alpha_fit_nu200_omega1"

omega_vect[] <- 1
fixed_parameter_list <- RSVsim_parameters(overrides = list("total_population" = 27536874,
                                                           "alpha_vect" = alpha_vect,
                                                           "nu" = 1 / 200,
                                                           "omega_vect" = omega_vect),
                                          contact_population_list = contact_population_list)

fixed_parameter_list <- apply_empirical_age_distribution(fixed_parameter_list)

##################
##### target #####
##################

# attack rate in 0 - 3 year olds
# annual incidence in adults aged 65+
# amplitude in everyone
# peak time (week) in everyone

annual_incidence_older_adults <- sum(cases_ti[ages_ti >= 65]) / pop_gr_65

target <- c(attack_rate[1:8], annual_incidence_older_adults, amplitude, peak_time)

###############################################################
##### summary function for calculating the model outputs ######
###############################################################

summary_fun <- local({
  pop_gr_65_local <- pop_gr_65
  
  function(out){
    
    out_subset <- out |> subset(age < 3 & time < 365.25)
    
    # cumulative incidence at 1, 2 and 3 years
    ar_y <- out_subset |>
      filter(age_chr %in% c("[0.92,1)", "[1.92,2)", "[2.92,3)")) |>
      group_by(age, age_chr, time) |>
      summarise(Sp = sum(Sp),
                total_pop = sum(Sp) + sum(Ep) + sum(Ip) + sum(Ss) + sum(Es) + sum(Is) + sum(R),
                .groups = "drop") |>
      dplyr::mutate(age_int = floor(age),
                    ps = 1 - Sp/total_pop) |> # ps is the proportion with at least 1 infection
      dplyr::group_by(age_int) |>
      dplyr::summarise(m_ps = mean(ps)) |>
      dplyr::arrange(age_int) |>
      dplyr::pull(m_ps)
    
    tot_inc <- out |> dplyr::group_by(week = floor(time / 7) + 1) |>
      dplyr::summarise(tot_incidence = sum(Incidence), .groups = "drop")
    
    # amplitude
    max_ <- max(max(tot_inc[,"tot_incidence"]), 1E-10)
    min_ <- max(min(tot_inc[,"tot_incidence"]), 1E-10)
    
    # annual incidence in adults aged 65+ approximated using total infections over one year
    inc_o <- sum(out[out$age >= 65 & out$time < 365.25, "Incidence"]) / pop_gr_65_local
    
    return(c(ar_y[1],
             ar_y[1],
             ar_y[2],
             ar_y[3],
             ar_y[1],
             ar_y[2],
             ar_y[1],
             ar_y[2],
             inc_o,
             (max_ - min_)/max_, # amplitude
             tot_inc$week[which.max(tot_inc$tot_incidence)][1] # peak time
    )
    )
    
  }
})

##########################
##### prior function #####
##########################

# uses latin hypercube sampling to sample from the prior space

prior_fun <- function(n_prior_attempts){
  
  x <- lhs::randomLHS(n_prior_attempts, 4)
  
  # adjusting the prior distributions
  x[,1] <- qunif(x[,1], min = 0.01, max = 0.2)
  x[,2] <- qunif(x[,2], min = 0.1, max = 0.3)
  x[,3] <- qunif(x[,3], min = 0, max = 365.25)
  x[,4] <- qunif(x[,4], min = 0.1, max = 0.5)
  
  return(as.matrix(x, nrow = n_prior_attempts))
}

##########################
##### prior function #####
##########################

prior_dens_fun <- function(x){
  return(c(dunif(x[1], min = 0.01, max = 0.2, log = FALSE),
           dunif(x[2], min = 0.1, max = 0.3, log = FALSE),
           dunif(x[3], min = 0, max = 365.25, log = FALSE),
           dunif(x[4], min = 0.1, max = 0.5, log = FALSE))
  )
}

##############################
##### distance function ######
##############################

dist_fun <- function(target, target_star){
  return(
    c(
      RSVsim_abs_dist_fun(head(target,-1), head(target_star, -1)),
      RSVsim_shortest_periodic_dist_fun(tail(target, 1), tail(target_star, 1), period = 365.25)
    )
  )
}

#################################
##### setting the tolerance #####
#################################
# calculating a tolerance that means at least 1 particle combination is accepted every 100 simulations
# getting 100 samples from the priors
set.seed(123)
n_check <- 200#1000
prior_params <- prior_fun(n_check)

# simulating the summary statistics for each particle
num_cores <- 9
cl <- makeCluster(num_cores)
registerDoParallel(cl)

clusterExport(cl, varlist = c("fixed_parameter_list", "time_step", "cohort_step_size", "summary_fun", "dist_fun", "target", "times_in", "warm_up_in", "RSVsim_run_model_fixed_age_distribution"))

prior_distances <- foreach(i = 1:n_check,
                           .packages = c("RSVsim", "dplyr", "tidyr")) %dopar% {
                             dist_fun(target,
                                      summary_fun(
                                        RSVsim_run_model_fixed_age_distribution(parameters = RSVsim_update_parameters(fixed_parameter_list, fitted_parameter_names, prior_params[i,]),
                                                                                times = times_in, # maximum time to run the model for
                                                                                cohort_step_size = cohort_step_size, # time at which to age people
                                                                                warm_up = warm_up_in)
                                      )
                             )
                           }
stopCluster(cl)

prior_distances <- do.call(rbind, prior_distances)

# create folder if it doesn't exist
if(!dir.exists(glue("outputs/{sn}"))){
  dir.create(glue("outputs/{sn}"))
}
saveRDS(prior_distances, file = glue("outputs/{sn}/prior_distances_100_min.rds"))

# calculating the number of particles for which all summary statistics are within the tolerance given different percentiles of the summary statistics
nsuccess <- rep(NA, n_check)
q <- seq(1/n_check, 1, 1/n_check)

for(i in 1:n_check){
  epsilon_check <- round(apply(prior_distances, 2, quantile, probs = c(q[i])), digits = 5)
  nsuccess[i] <- sum(sapply(1:n_check, function(j){all(prior_distances[j,] <= epsilon_check)}))
}

acceptance_rate <- sort(c(0.1, 0.5, 1, 5, 10, seq(20, 70, 25)), decreasing = TRUE)

q_percentile <- vector(mode = "list", length = length(acceptance_rate))

q_percentile <- lapply(acceptance_rate, function(ar){q[min(which(nsuccess/n_check * 100 > ar))]}) # 1

q_percentile_ar1 <- q_percentile[[which(acceptance_rate == 1)]]

# selecting a tolerance where at least 1 particle is accepted for the 100 simulations
epsilon_matrix <- as.matrix(do.call(rbind, lapply(q_percentile, function(q){round(apply(prior_distances, 2, quantile, probs = q), digits = 3)})))

epsilon <- epsilon_matrix[which(acceptance_rate == 1),]

G <- nrow(epsilon_matrix)

###############################################
##### running the ABC-rejection algorithm #####
###############################################

nparticles = 200

used_seeds_all <- seq(1, nparticles)

ncores <- 9

# fit <- RSVsim_ABC_rejection(target = target,
#                             epsilon = epsilon,
#                             summary_fun = summary_fun,
#                             dist_fun = dist_fun,
#                             prior_fun = prior_fun,
#                             n_prior_attempts = 10000,
#                             nparticles = nparticles,
#                             used_seeds_all = used_seeds_all,
#                             ncores = ncores,
#                             fitted_parameter_names = fitted_parameter_names,
#                             fixed_parameter_list = fixed_parameter_list,
#                             times = seq(0, 365.25*11 + 5.75 - time_step, time_step), # maximum time to run the model for
#                             cohort_step_size = cohort_step_size, # time at which to age people\
#                             warm_up = 10*365.25)
#
# saveRDS(fit, file = "model_fit_ABC_reject.rds")

#########################################
##### running the ABC-SMC algorithm #####
#########################################

used_seed_matrix <- matrix(seq(1, nparticles * G), nrow = G)

fit_smc <- RSVsim_ABC_SMC_fixed_age_distribution(target = target,
                                                 epsilon_matrix = epsilon_matrix,
                                                 summary_fun = summary_fun,
                                                 dist_fun = dist_fun,
                                                 prior_fun = prior_fun,
                                                 n_prior_attempts = n_prior_attempts,
                                                 used_seed_matrix = used_seed_matrix,
                                                 prior_dens_fun = prior_dens_fun,
                                                 particle_low = c(0, 0.1, 0, 0.1),
                                                 particle_up = c(0.2, 0.3, 365.25, 0.5),
                                                 nparticles = nparticles,
                                                 ncores = ncores,
                                                 fitted_parameter_names = fitted_parameter_names,
                                                 fixed_parameter_list = fixed_parameter_list,
                                                 times = times_in, # maximum time to run the model for
                                                 cohort_step_size = cohort_step_size, # time at which to age people
                                                 warm_up = warm_up_in)

saveRDS(fit_smc, file = glue("outputs/{sn}/model_fit_ABC_SMC_100_min.rds"))
