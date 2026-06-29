rm(list = ls())
source("data/data.r")
source("R/set_up_variables.r")

sn <- "alpha_fit_nu200_omega1"
omega_vect[] <- 1
fixed_parameter_list <- RSVsim_parameters(overrides = list("total_population" = 27536874,
                                                           "alpha_vect" = alpha_vect,
                                                           "nu" = 1 / 200,
                                                           "omega_vect" = omega_vect),
                                          contact_population_list = contact_population_list)
fixed_parameter_list <- apply_empirical_age_distribution(fixed_parameter_list)

fit_all <- readRDS(file = glue("outputs/{sn}/model_fit_ABC_SMC_100_min.rds"))
G <- length(fit_all$fitted_parameters)

fit <- fit_all$fitted_parameters[[G]]
weights <- fit_all$weights[[G]]

colnames(fit) <- c("b0", "b1", "phi", "alpha")
fit <- as.data.frame(fit)
fit$particle_number <- 1:nrow(fit)

# posterior parameters
fit[,"mean_b0"] <- rep(NA, nrow(fit))
fit[,"mean_b1"] <- rep(NA, nrow(fit))
fit[,"mean_phi"] <- rep(NA, nrow(fit))
fit[,"mean_alpha"] <- rep(NA, nrow(fit))

for(i in 1:nrow(fit)){
  fit[i, "mean_b0"] <- mean(fit[1:i, "b0"])
  fit[i, "mean_b1"] <- mean(fit[1:i, "b1"])
  fit[i, "mean_phi"] <- mean(fit[1:i, "phi"])
  fit[i, "mean_alpha"] <- mean(fit[1:i, "alpha"])
}

plot(fit$particle_number, fit$mean_b0, type = "l")
plot(fit$particle_number, fit$mean_b1, type = "l")
plot(fit$particle_number, fit$mean_phi, type = "l")
plot(fit$particle_number, fit$mean_alpha, type = "l")

plot(fit$particle_number, fit$b0, type = "l")
plot(fit$particle_number, fit$b1, type = "l")
plot(fit$particle_number, fit$phi, type = "l")
plot(fit$particle_number, fit$alpha, type = "l")

######################
##### model runs #####
######################

num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, varlist = c("fixed_parameter_list", "time_step", "fit", "cohort_step_size", "times_in", "warm_up_in"))
clusterExport(cl, varlist = c("RSVsim_run_model_fixed_age_distribution"))

model_runs <- foreach(i = 1:nrow(fit),
                      .packages = c("RSVsim", "dplyr", "tidyr")) %dopar% {
                        
                        RSVsim_run_model_fixed_age_distribution(parameters = RSVsim_update_parameters(fixed_parameter_list,
                                                                                                      fitted_parameter_names,
                                                                                                      c(fit[i, "b0"], fit[i, "b1"], fit[i, "phi"], fit[i, "alpha"])),
                                                                 times = times_in, # maximum time to run the model for
                                                                 cohort_step_size = cohort_step_size, # time at which to age people
                                                                 warm_up = warm_up_in) |> mutate(prior_function_seed = i) |>
                          filter(vacc_state == 1) |>
                          select(-VaccState_rate, -doses)
                      }
stopCluster(cl)

saveRDS(model_runs, file = glue("outputs/{sn}/model_runs_100_min.rds"))
model_runs <- readRDS(glue("outputs/{sn}/model_runs_100_min.rds"))

#####################
##### incidence #####
#####################

inc_df <- lapply(model_runs, function(df){
  df |> dplyr::mutate(week = floor(time / 7) + 1) |> group_by(week) |>
    summarise(tot_inc = sum(Incidence))  |>
    mutate(sc_inc = tot_inc / max(tot_inc)*1.2)}) |> bind_rows() |>
  group_by(week) |>
  summarise(median_tot = median(tot_inc),
            lower_tot = quantile(tot_inc, probs = lower),
            upper_tot = quantile(tot_inc, probs = upper),
            median_sc = median(sc_inc),
            lower_sc = quantile(sc_inc, probs = lower),
            upper_sc = quantile(sc_inc, probs = upper))


age_groups_df <- data.frame(age_group_plot = floor(contact_population_list$age_limits / 5) * 5,
                            pop = contact_population_list$age_distribution * fixed_parameter_list$total_population) |>
  group_by(age_group_plot) |> summarise(pop = sum(pop))


inc_by_age_df <- lapply(model_runs,
                        function(df){
                          suppressMessages(df |> mutate(week = floor(time / 7) + 1,
                                                        age_group_plot = floor(age / 5) * 5) |>
                                             group_by(age_group_plot, week, prior_function_seed) |>
                                             summarise(,
                                                       Incidence = sum(Incidence),
                                                       DetIncidence = sum(DetIncidence),
                                             ) |> ungroup() |> left_join(age_groups_df, by = "age_group_plot")
                          )
                        }
)

inc_by_age_df <- inc_by_age_df |> bind_rows() |> mutate(age_group_plot = ifelse(age_group_plot != 75, paste0("[",age_group_plot,",",age_group_plot + 5,")"), "[75, 90)"))

inc_by_age_df$age_group_plot <- factor(inc_by_age_df$age_group_plot, levels = c(paste0("[",seq(0, 70, 5), ",",seq(5, 75, 5), ")"), "[75, 90)"))

inc_by_age_df_sum <- inc_by_age_df |> mutate(Incidence = Incidence / pop * 1000)  |> group_by(age_group_plot, week) |>
  summarise(mI = median(Incidence),
            lI = quantile(Incidence, probs = lower),
            uI = quantile(Incidence, probs = upper))

max_week_df <- inc_by_age_df |> group_by(age_group_plot, prior_function_seed) |> filter(Incidence == max(Incidence)) |> ungroup() |> group_by(age_group_plot) |>
  summarise(mW = median(week), lW = quantile(week, probs = lower), uW = quantile(week, probs = upper))

saveRDS(inc_df, glue("outputs/{sn}/inc_df.rds"))
saveRDS(inc_by_age_df, glue("outputs/{sn}/inc_by_age_df.rds"))
saveRDS(inc_by_age_df_sum,glue( "outputs/{sn}/inc_by_age_df_sum.rds"))
saveRDS(max_week_df, glue("outputs/{sn}/max_week_df.rds"))

########################
##### attack rates #####
########################

# in older adults sum of incidence / number of people in that age group
ar_df <- inc_by_age_df |> group_by(age_group_plot, prior_function_seed, pop) |> summarise(AR = sum(Incidence)) |> mutate(AR = AR / pop) |> ungroup() |>
  group_by(age_group_plot) |>
  summarise(mAR = median(AR),
            lAR = quantile(AR, probs = lower),
            uAR = quantile(AR, probs = upper))

out_subset <- lapply(model_runs, function(df){
  df |> subset(age < 3 & time < 365.25) |>
    arrange(age, time) |>
    mutate(cohort_birth = ifelse(
      round(time, digits = 7) %in% age_times,
      round(age - time / 365.25, digits = 7),
      NA)) |>
    fill(cohort_birth, .direction = "down") |> na.omit() |>
    mutate(age_continuous = cohort_birth + time / 365.25) |>
    subset(age_continuous > 0) |>
    mutate(ar = 1 - (Sp/(Sp + Ep + Ip + Ss + Es + Is + R))) |>
    select(age, time, ar, cohort_birth, age_continuous, prior_function_seed)
}
) |> bind_rows()


ar_df_03 <- out_subset %>%
  filter(age %in% c(0.91667, 1.91667, 2.91667)) %>%
  group_by(age, prior_function_seed) %>%
  summarise(ar = mean(ar)) %>%
  ungroup() %>%
  group_by(age) %>%
  summarise(ar_median = median(ar),
         ar_upper = quantile(ar, 0.975),
         ar_lower = quantile(ar, 0.025))

ar_df_03

ar_df_c <- out_subset |>
  group_by(cohort_birth, age_continuous) |>
  summarise(ar_m = median(ar),
            ar_l = quantile(ar, lower),
            ar_u = quantile(ar, upper))

saveRDS(ar_df, glue("outputs/{sn}/ar_df.rds"))
saveRDS(ar_df_c, glue("outputs/{sn}/ar_df_c.rds"))
saveRDS(ar_df_03, glue("outputs/{sn}/ar_df_03.rds"))
