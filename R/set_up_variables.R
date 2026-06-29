library(RSVsim)
library(tidyverse)
source("R/run_model_fixed_age_distribution.R")
source("R/ABC_SMC_fixed_age_distribution.R")

############################
##### set up variables #####
############################

step <- 1/12
cohort_step_size <- step * 365.25

age_times <- seq(0, 365.25 - cohort_step_size, cohort_step_size)

age_limits <- c(seq(0, 5, step), seq(10, 75, 5))

contact_population_list <- RSVsim_contact_matrix(country = "United Kingdom", age_limits = age_limits)

time_step <- 0.25

times_in <- seq(0, 365.25*11 + 5.75 - time_step, time_step)
#times_in <- seq(0, 365.25* +20+ 5.75 - time_step, time_step)

warm_up_in <- 365.25 * 10
#warm_up_in <- 365.25*10
n_prior_attempts <- 10000

nAges <- length(age_limits)

# fit a single alpha value across all ages via alpha_vect[]
alpha_vect <- rep(0.3, nAges)

omega_vect <- rep(NA, nAges)

omega_vect[] <- 1

fixed_parameter_list <- RSVsim_parameters(overrides = list("total_population" = 27536874,
                                                           "alpha_vect" = alpha_vect,
                                                           "nu" = 1 / 200,
                                                           "omega_vect" = omega_vect),
                                          contact_population_list = contact_population_list)

# On the fixed-demography branch we initialise the model with the empirical
# age distribution from the contact matrix object. The custom model runner then
# applies an explicit constant demographic flow so these age-class totals remain
# effectively fixed over time.
apply_empirical_age_distribution <- function(parameter_list, infected_prop = 0.001) {
  age_distribution_population <- contact_population_list$age_distribution * parameter_list$total_population
  parameter_list$Sp0[,] <- 0
  parameter_list$Ep0[,] <- 0
  parameter_list$Ip0[,] <- 0
  parameter_list$Ss0[,] <- 0
  parameter_list$Es0[,] <- 0
  parameter_list$Is0[,] <- 0
  parameter_list$R0[,] <- 0
  parameter_list$Sp0[, 1] <- age_distribution_population * (1 - infected_prop)
  parameter_list$Ip0[, 1] <- age_distribution_population * infected_prop
  parameter_list
}

fixed_parameter_list <- apply_empirical_age_distribution(fixed_parameter_list)

fitted_parameter_names <- c("b0", "b1", "phi", "alpha_vect[]")

alpha_lower <- 0.1
alpha_upper <- 0.5

lower <- 0.025
upper <- 0.975

# population sizes
pop_all <- contact_population_list$age_distribution * fixed_parameter_list$total_population
pop_gr_65 <- sum(pop_all[which(age_limits >= 65)])
pop_gr_70 <- sum(pop_all[which(age_limits >= 70)])

# age specific number of hospitalisations per person
hosp_rate <- c(5338.91375, 1751.145228, 715.9402014, 208.4621667, 10.86702733, 2.58129891, 2.321354938, 1.868628613,
               1.320440418, 1.834551244, 1.883359955, 2.861263391, 4.744902927, 6.645799647, 9.348531561, 14.3538506,
               22.62406397, 32.94296492, 87.63366855) / 100000

age_limits_hosp_rate <- c(0, 0.5, 1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75) # lower age limits of the hospitalisation rate data
