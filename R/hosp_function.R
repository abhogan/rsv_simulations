RSVsim_hospitalisations <- function(prob_hosp_given_inf,
                                    age_limits_hosp_rate,
                                    sim_scenario,
                                    VE_hosp_given_inf,
                                    VE_hosp_vacc_states){
  
  age_limits <- unique(sim_scenario$age)
  
  age_limits <- round(age_limits, digits = 5)
  
  age_limits_hosp_rate <- round(age_limits_hosp_rate, digits = 5)
  
  sim_scenario[,"time"] <- round(sim_scenario[, "time"], digits = 5)
  
  if(base::is.unsorted(age_limits_hosp_rate)){
    stop("RSVsim_hospitalisations: age_limits_hosp_rate must be in ascending order and hosp_rate values must be in the same position as the corresponding age")
  }
  
  if(!all(age_limits_hosp_rate %in% age_limits)){
    stop("RSVsim_hospitalisations: not all age_limits_hosp_rate are present in sim_scenario")
  }
  
  if(nrow(VE_hosp_given_inf) != length(age_limits) || ncol(VE_hosp_given_inf) != length(VE_hosp_vacc_states)){
    stop("RSVsim_hospitalisations: VE_hosp_given_inf must be a matrix with the rows corresponding to the ages in age_limits
         and the columns corresponding to the VE_hosp_vacc_states")
  }
  
  if(length(prob_hosp_given_inf) != length(age_limits_hosp_rate)){
    stop("RSVsim_hospitalisations: prob_hosp_given_inf must be the same length as age_limits_hosp_rate")
  }
  
  age_limits_index <- data.table::data.table(age_index = base::findInterval(age_limits, age_limits_hosp_rate),
                                             age = age_limits)
  age_limits_index[, age_hosp := age_limits_hosp_rate[age_index]]
  
  sim_scenario <- data.table::as.data.table(sim_scenario)
  sim_scenario[age_limits_index, on = .(age), age_hosp := i.age_hosp]
  
  inc_per_person_df <- data.table::data.table(prob_hosp_given_inf_uv = prob_hosp_given_inf,
                                              age_hosp = age_limits_hosp_rate)
  
  VE_hosp_df <- data.table::as.data.table(VE_hosp_given_inf)
  data.table::setnames(VE_hosp_df, as.character(VE_hosp_vacc_states))
  VE_hosp_df[ , age := age_limits]
  VE_hosp_df <- data.table::melt(VE_hosp_df,
                                 id.vars = "age",
                                 variable.name = "vacc_state",
                                 value.name = "VE_hosp_given_inf")
  VE_hosp_df[, vacc_state := as.numeric(as.character(vacc_state))]
  
  sim_scenario[inc_per_person_df, on = .(age_hosp), prob_hosp_given_inf_uv := i.prob_hosp_given_inf_uv]
  sim_scenario[VE_hosp_df, on = .(vacc_state, age), VE_hosp_given_inf := i.VE_hosp_given_inf]
  sim_scenario[, ":="(Hospitalisations = Incidence * (1 - data.table::fcoalesce(VE_hosp_given_inf, 0)) * prob_hosp_given_inf_uv,
                      prob_hosp_given_inf_uv = NULL,
                      VE_hosp_given_inf = NULL,
                      age_hosp = NULL)]
  
  return(as.data.frame(sim_scenario))
}
