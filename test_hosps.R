
times <- seq(0, 365.25*11, 0.25)
warm_up <- 365.25 * 10
cohort_step_size <- 1/12 * 365.25

age_limits <- c(seq(0, 5, 1/12), seq(10, 75, 5))

contact_population_list <- RSVsim_contact_matrix(age_limits = age_limits)

parameters <- RSVsim_parameters(contact_population_list = contact_population_list)

sim_status_quo <- RSVsim_run_model(parameters = parameters,
                                   times = times,
                                   cohort_step_size = cohort_step_size,
                                   warm_up = warm_up)



# vaccine efficacy against hospitalisation conditional on infection
VE_hosp <- matrix(rep(0.8, length(age_limits)), ncol = 1) # age-specific vaccine efficacy against hospitalisation
VE_inf <- parameters$VE_inf
VE_hosp_given_inf <- RSVsim_VE_hosp_given_inf(VE_hosp = VE_hosp, VE_inf = VE_inf)

prob_hosp_given_inf_list <- RSVsim_prob_hosp_given_inf(hosp_rate = hosp_rate,
                                                       age_limits_hosp_rate = age_limits_hosp_rate,
                                                       hosp_rate_min_time = 0, # time frame over which to sum incidence
                                                       hosp_rate_max_time = 365.25,
                                                       sim_status_quo = sim_status_quo,
                                                       VE_hosp_given_inf = VE_hosp_given_inf,
                                                       VE_hosp_vacc_states = c(2) # vector identifying the vaccination states that are vaccinated (position corresponds to the column position in VE_hosp)
)
