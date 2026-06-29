RSVsim_run_model_fixed_age_distribution <- function(parameters,
                                                    times = seq(0, 365.25 * 1, 0.25),
                                                    cohort_step_size = 1 / 12 * 365.25,
                                                    warm_up = NULL) {
  # This branch uses an explicit stationary demographic flow rather than the
  # package default ageing step. A constant number of people enters the first
  # age class each cohort step, the same number progresses across every age
  # boundary proportionally across compartments, and the same number exits the
  # final age class. This preserves the empirical age distribution and total
  # population size over time.
  max_t <- max(times)
  RSV_dust <- invisible(dust2::dust_system_create(generator = getFromNamespace("RSV_ODE", "RSVsim"),
                                                  pars = parameters,
                                                  n_particles = 1,
                                                  n_groups = 1,
                                                  time = 0,
                                                  deterministic = TRUE))
  dust2::dust_system_set_state_initial(RSV_dust)
  states <- dust2::dust_unpack_index(RSV_dust)
  states_order <- names(states)
  pop_states <- c("Sp", "Ep", "Ip", "Ss", "Es", "Is", "R")
  target_age_totals <- rowSums(parameters$Sp0 + parameters$Ep0 + parameters$Ip0 +
                                 parameters$Ss0 + parameters$Es0 + parameters$Is0 + parameters$R0)
  states <- as.character(rep(states_order, lengths(states)))
  n_states <- length(states_order)
  ages <- base::round(rep(rep(parameters$age_limits, parameters$nVaccStates), n_states), digits = 5)
  vacc_states <- base::round(rep(rep(1:parameters$nVaccStates, each = parameters$nAges), n_states), digits = 5)
  age_chr <- rep(rep(parameters$age_chr, parameters$nVaccStates), n_states)
  output_variable_names <- states_order[!(states_order %in% c("Sp", "Ep", "Ip", "Ss", "Es", "Is", "R"))]
  dust_df <- data.frame(state = as.character(states),
                        age = ages,
                        age_chr = age_chr,
                        vacc_state = vacc_states,
                        index = 1:length(states),
                        stringsAsFactors = FALSE)
  
  if (is.numeric(cohort_step_size)) {
    birth_flow <- target_age_totals[1] * cohort_step_size / parameters$size_cohorts[1]
    transition_rate <- birth_flow / target_age_totals
    if (base::round(max(diff(times)), digits = 5) >= base::round(cohort_step_size, digits = 5)) {
      stop("The maximum time difference is greater than or equal to the cohort step size")
    }
    if (any(transition_rate > 1 + 1e-8)) {
      stop("RSVsim_run_model_fixed_age_distribution: cohort_step_size is too large for the explicit demographic flow")
    }
    
    dust_df <- as.data.frame(dplyr::ungroup(
      dplyr::mutate(
        dplyr::group_by(dplyr::arrange(dust_df, index), state, vacc_state),
        cohort_ageing = dplyr::if_else(state %in% output_variable_names, 0, 1),
        lag_index = dplyr::if_else(cohort_ageing == 0, index, dplyr::lag(index, 1))
      )
    ))
    
    age_index <- match(dust_df$age, base::round(parameters$age_limits, digits = 5))
    transition_rate_all <- transition_rate[age_index]
    births_zero_index <- which(is.na(dust_df$lag_index))
    birth_index <- which(is.na(dust_df$lag_index) &
                           dust_df$age == 0 &
                           dust_df$state == "Sp" &
                           dust_df$vacc_state == 1)
    
    n_steps <- ceiling(max_t / cohort_step_size)
    out_list <- vector(mode = "list", length = n_steps)
    times_all <- sort(unique(base::round(c(times, 1:n_steps * cohort_step_size), digits = 5)))
    times_in <- lapply(1:n_steps, FUN = function(i) {
      c(times_all[times_all >= base::round(((i - 1) * cohort_step_size), digits = 5) &
                    times_all < base::round((i * cohort_step_size), digits = 5)])
    })
    if (max(times_in[[n_steps]]) < max_t) {
      times_in[[n_steps]] <- c(times_in[[n_steps]], max_t)
    }
    
    for (i in 1:n_steps) {
      out <- t(dust2::dust_system_simulate(RSV_dust, times = times_in[[i]]))
      next_state <- out[base::length(times_in[[i]]), ]
      out <- base::cbind(out, times_in[[i]])
      out_list[[i]] <- out
      
      transition_ct <- next_state * transition_rate_all
      lag_transition_ct <- transition_ct[dust_df$lag_index]
      lag_transition_ct[births_zero_index] <- 0
      next_value <- next_state - transition_ct + lag_transition_ct
      next_value[birth_index] <- next_value[birth_index] + birth_flow
      
      dust2::dust_system_set_state(sys = RSV_dust, state = next_value)
    }
    out_checkout <- as.data.frame(do.call(rbind, out_list))
    rm(out_list)
  } else {
    times <- base::round(times, digits = 5)
    out_checkout <- as.data.frame(t(dust2::dust_system_simulate(RSV_dust, times = times)))
    out_checkout[, "time"] <- times
  }
  
  colnames(out_checkout) <- c(base::as.character(dust_df[, "index"]), "time")
  if (!is.null(warm_up)) {
    out_checkout <- dplyr::filter(out_checkout, time >= times[min(which(times >= warm_up)) - 1])
  }
  data.table::setDTthreads(1)
  out_checkout <- data.table::as.data.table(out_checkout)
  dust_df <- data.table::as.data.table(dust_df)
  out_checkout <- data.table::melt(out_checkout, id.vars = "time", variable.name = "index", value.name = "value")
  out_checkout[, `:=`(index, as.integer(as.character(index)))]
  out_checkout <- out_checkout[dust_df[, c("state", "age", "age_chr", "vacc_state", "index")], on = "index", nomatch = NA]
  out_checkout[, `:=`(age = round(age, 5), vacc_state = as.integer(vacc_state), time = round(time, 5))]
  out_checkout <- data.table::dcast(out_checkout, age + age_chr + vacc_state + time ~ state, value.var = "value", drop = TRUE)
  
  pop_check <- base::as.data.frame(out_checkout[, .(total = sum(rowSums(.SD))), by = time, .SDcols = c("Sp", "Ep", "Ip", "Ss", "Es", "Is", "R")])
  data.table::setorder(out_checkout, vacc_state, age, time)
  data.table::setcolorder(out_checkout, c("age", "age_chr", "vacc_state", "time", states_order))
  out_checkout <- dplyr::ungroup(
    dplyr::mutate(
      dplyr::group_by(as.data.frame(out_checkout), age_chr, vacc_state),
      Incidence = tidyr::replace_na(Incidence - dplyr::lag(Incidence, 1), 0),
      DetIncidence = tidyr::replace_na(DetIncidence - dplyr::lag(DetIncidence, 1), 0),
      doses = tidyr::replace_na(doses - dplyr::lag(doses, 1), 0)
    )
  )
  if (!is.null(warm_up)) {
    out_checkout <- dplyr::mutate(dplyr::filter(out_checkout, time >= warm_up), time = time - warm_up)
  }
  max_pop_dev <- max(abs(pop_check[, "total"] - parameters$total_population))
  if (max_pop_dev > 0.1) {
    stop(sprintf("RSVsim_run_model_fixed_age_distribution: population does not sum to the correct number (max deviation = %.6f)", max_pop_dev))
  }
  return(as.data.frame(out_checkout))
}
