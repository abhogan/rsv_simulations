rm(list = ls())
source("data/data.r")
source("R/set_up_variables.r")
source("R/hosp_function.R")

# Run baseline and vaccine scenarios using the updated fitted parameters from the
# uniform-age calibration workflow on this branch.
sn <- "alpha_fit_nu200_omega1"

omega_vect[] <- 1

dir.create(glue("outputs/{sn}"), recursive = TRUE, showWarnings = FALSE)

fixed_parameter_list <- RSVsim_parameters(
  overrides = list(
    total_population = 27536874,
    nu = 1 / 200,
    omega_vect = omega_vect
  ),
  contact_population_list = contact_population_list
)

fit_all <- readRDS(file = glue("outputs/{sn}/model_fit_ABC_SMC_100_min.rds"))
G <- length(fit_all$fitted_parameters)
fit <- fit_all$fitted_parameters[[G]] |>
  as.data.frame()

fit <- fit[1:20,]

colnames(fit) <- c("b0", "b1", "phi", "alpha")
fit$prior_function_seed <- seq_len(nrow(fit))

VE_hosp <- matrix(rep(0.8, length(age_limits)), ncol = 1)

run_with_fit <- function(parameter_list, particle_row) {
  RSVsim_run_model(
    parameters = RSVsim_update_parameters(
      parameter_list,
      fitted_parameter_names,
      c(particle_row$b0, particle_row$b1, particle_row$phi, particle_row$alpha)
    ),
    times = times_in,
    cohort_step_size = cohort_step_size,
    warm_up = warm_up_in
  )
}

get_prob_hosp_given_inf <- function(sim_status_quo, parameter_list) {
  VE_hosp_given_inf <- RSVsim_VE_hosp_given_inf(
    VE_hosp = VE_hosp,
    VE_inf = parameter_list$VE_inf
  )

  RSVsim_prob_hosp_given_inf(
    hosp_rate = hosp_rate,
    age_limits_hosp_rate = age_limits_hosp_rate,
    hosp_rate_min_time = 0,
    hosp_rate_max_time = 365.25,
    sim_status_quo = sim_status_quo,
    VE_hosp_given_inf = VE_hosp_given_inf,
    VE_hosp_vacc_states = c(2)
  )$prob_hosp_given_inf_uv
}

add_hospitalisations <- function(sim_scenario, prob_hosp_given_inf, parameter_list) {
  VE_hosp_given_inf <- RSVsim_VE_hosp_given_inf(
    VE_hosp = VE_hosp,
    VE_inf = parameter_list$VE_inf
  )

  RSVsim_hospitalisations(
    prob_hosp_given_inf = prob_hosp_given_inf,
    age_limits_hosp_rate = age_limits_hosp_rate,
    sim_scenario = sim_scenario,
    VE_hosp_given_inf = VE_hosp_given_inf,
    VE_hosp_vacc_states = c(2)
  )
}

summarise_totals <- function(model_runs, scenario_name) {
  lapply(model_runs, function(df) {
    df |>
      mutate(week = floor(time / 7) + 1) |>
      group_by(week) |>
      summarise(
        tot_inc = sum(Incidence),
        tot_hosp = sum(Hospitalisations),
        .groups = "drop"
      ) |>
      mutate(sc_inc = ifelse(max(tot_inc) > 0, tot_inc / max(tot_inc), 0))
  }) |>
    bind_rows() |>
    group_by(week) |>
    summarise(
      median_tot = median(tot_inc),
      lower_tot = quantile(tot_inc, probs = lower),
      upper_tot = quantile(tot_inc, probs = upper),
      median_sc = median(sc_inc),
      lower_sc = quantile(sc_inc, probs = lower),
      upper_sc = quantile(sc_inc, probs = upper),
      median_hosp = median(tot_hosp),
      lower_hosp = quantile(tot_hosp, probs = lower),
      upper_hosp = quantile(tot_hosp, probs = upper),
      .groups = "drop"
    ) |>
    mutate(scenario = scenario_name)
}

summarise_by_age <- function(model_runs, age_groups_df, scenario_name) {
  out <- lapply(model_runs, function(df) {
    suppressMessages(
      df |>
        mutate(
          week = floor(time / 7) + 1,
          age_group_plot = floor(age / 5) * 5
        ) |>
        group_by(age_group_plot, week, prior_function_seed) |>
        summarise(
          Incidence = sum(Incidence),
          DetIncidence = sum(DetIncidence),
          Doses = sum(doses),
          Hospitalisations = sum(Hospitalisations),
          .groups = "drop"
        ) |>
        left_join(age_groups_df, by = "age_group_plot") |>
        mutate(scenario = scenario_name)
    )
  }) |>
    bind_rows() |>
    mutate(age_group_plot = ifelse(age_group_plot != 75,
                                   paste0("[", age_group_plot, ",", age_group_plot + 5, ")"),
                                   "[75, 90)"))

  out$age_group_plot <- factor(
    out$age_group_plot,
    levels = c(paste0("[", seq(0, 70, 5), ",", seq(5, 75, 5), ")"), "[75, 90)")
  )

  out
}

################################
##### model runs: find IHR #####
################################

num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(
  cl,
  varlist = c(
    "VE_hosp",
    "age_limits_hosp_rate",
    "cohort_step_size",
    "fit",
    "fitted_parameter_names",
    "fixed_parameter_list",
    "get_prob_hosp_given_inf",
    "hosp_rate",
    "run_with_fit",
    "times_in",
    "warm_up_in"
  )
)

prob_hosps <- foreach(i = seq_len(nrow(fit)), .packages = c("RSVsim", "dplyr", "tidyr")) %dopar% {
  sim_status_quo <- run_with_fit(fixed_parameter_list, fit[i, ])
  get_prob_hosp_given_inf(sim_status_quo, fixed_parameter_list)
}
stopCluster(cl)

################################
##### model runs: baseline #####
################################

num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(
  cl,
  varlist = c(
    "add_hospitalisations",
    "age_limits_hosp_rate",
    "cohort_step_size",
    "fit",
    "fixed_parameter_list",
    "prob_hosps",
    "run_with_fit",
    "times_in",
    "VE_hosp",
    "warm_up_in"
  )
)

model_runs_baseline <- foreach(i = seq_len(nrow(fit)), .packages = c("RSVsim", "dplyr", "tidyr")) %dopar% {
  sim_baseline <- run_with_fit(fixed_parameter_list, fit[i, ]) |>
    mutate(prior_function_seed = i)

  add_hospitalisations(sim_baseline, prob_hosps[[i]], fixed_parameter_list)
}
stopCluster(cl)

saveRDS(model_runs_baseline, file = glue("outputs/{sn}/model_runs_baseline.rds"))

################################
##### model runs: vaccine ######
################################

vaccine_times <- c(0, c(60, 180) + warm_up_in)
vaccine_period <- diff(c(vaccine_times, max(times_in)))
nVaccTimes <- length(vaccine_times)

vaccine_cov <- rbind(
  matrix(rep(0, (nAges - 3) * nVaccTimes), nrow = (nAges - 3), ncol = nVaccTimes),
  c(0, 0.80, 0),
  c(0, 0.80, 0),
  c(0, 0.80, 0)
)

nVaccStates <- 2
VE_inf <- matrix(rep(0.5, nAges * (nVaccStates - 1)), nrow = nAges)
gamma_vaccine <- 1 / (365.25 * 2)

parameters_vaccine <- RSVsim_parameters(
  overrides = list(
    vaccine_times = vaccine_times,
    vaccine_period = vaccine_period,
    nVaccTimes = nVaccTimes,
    vaccine_cov = vaccine_cov,
    nVaccStates = nVaccStates,
    VE_inf = VE_inf,
    gamma_vaccine = gamma_vaccine,
    total_population = 27536874,
    alpha_vect = alpha_vect,
    nu = 1 / 200,
    omega_vect = omega_vect
  ),
  contact_population_list = contact_population_list
)

num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(
  cl,
  varlist = c(
    "add_hospitalisations",
    "cohort_step_size",
    "fit",
    "parameters_vaccine",
    "prob_hosps",
    "run_with_fit",
    "times_in",
    "VE_hosp",
    "warm_up_in"
  )
)

model_runs_vaccine <- foreach(i = seq_len(nrow(fit)), .packages = c("RSVsim", "dplyr", "tidyr")) %dopar% {
  sim_vaccine <- run_with_fit(parameters_vaccine, fit[i, ]) |>
    mutate(prior_function_seed = i)

  add_hospitalisations(sim_vaccine, prob_hosps[[i]], parameters_vaccine)
}
stopCluster(cl)

saveRDS(model_runs_vaccine, file = glue("outputs/{sn}/model_runs_vaccine.rds"))

###############################################################################

inc_df_baseline <- summarise_totals(model_runs_baseline, "baseline")
inc_df_vaccine <- summarise_totals(model_runs_vaccine, "vaccine")

age_groups_df <- data.frame(
  age_group_plot = floor(contact_population_list$age_limits / 5) * 5,
  pop = contact_population_list$age_distribution * fixed_parameter_list$total_population
) |>
  group_by(age_group_plot) |>
  summarise(pop = sum(pop), .groups = "drop")

inc_by_age_df_baseline <- summarise_by_age(model_runs_baseline, age_groups_df, "baseline")
inc_by_age_df_vaccine <- summarise_by_age(model_runs_vaccine, age_groups_df, "vaccine")

inc_by_age_df_baseline_sum <- inc_by_age_df_baseline |>
  mutate(Incidence = Incidence / pop * 1000) |>
  group_by(age_group_plot, week, scenario) |>
  summarise(
    mI = median(Incidence),
    lI = quantile(Incidence, probs = lower),
    uI = quantile(Incidence, probs = upper),
    mH = median(Hospitalisations),
    lH = quantile(Hospitalisations, probs = lower),
    uH = quantile(Hospitalisations, probs = upper),
    .groups = "drop"
  )

inc_by_age_df_vaccine_sum <- inc_by_age_df_vaccine |>
  mutate(Incidence = Incidence / pop * 1000) |>
  group_by(age_group_plot, week, scenario) |>
  summarise(
    mI = median(Incidence),
    lI = quantile(Incidence, probs = lower),
    uI = quantile(Incidence, probs = upper),
    mH = median(Hospitalisations),
    lH = quantile(Hospitalisations, probs = lower),
    uH = quantile(Hospitalisations, probs = upper),
    .groups = "drop"
  )

inc_df <- bind_rows(inc_df_baseline, inc_df_vaccine)
inc_by_age_df <- bind_rows(inc_by_age_df_baseline, inc_by_age_df_vaccine)
inc_by_age_df_sum <- bind_rows(inc_by_age_df_baseline_sum, inc_by_age_df_vaccine_sum)

saveRDS(inc_df, glue("outputs/{sn}/inc_df_vaccine.rds"))
saveRDS(inc_by_age_df, glue("outputs/{sn}/inc_by_age_df_vaccine.rds"))
saveRDS(inc_by_age_df_sum, glue("outputs/{sn}/inc_by_age_df_sum_vaccine.rds"))
