rm(list = ls())
library(tidyverse)
library(glue)

sn <- "alpha_fit_nu200_omega1"
total_pop <- 27536874

dir.create(glue("outputs/{sn}/figures"), recursive = TRUE, showWarnings = FALSE)

inc_df <- readRDS(glue("outputs/{sn}/inc_df_vaccine.rds"))
inc_by_age_df <- readRDS(glue("outputs/{sn}/inc_by_age_df_vaccine.rds"))
inc_by_age_df_sum <- readRDS(glue("outputs/{sn}/inc_by_age_df_sum_vaccine.rds"))

hosp_time_plot <- ggplot(data = inc_df, aes(x = week, y = median_hosp / total_pop * 100000, colour = scenario)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = lower_hosp / total_pop * 100000,
      ymax = upper_hosp / total_pop * 100000,
      fill = scenario
    ),
    alpha = 0.15,
    colour = NA
  ) +
  theme_minimal() +
  labs(x = "Week", y = "Hospitalisations per 100,000") +
  scale_colour_manual(values = c(baseline = "red3", vaccine = "bisque3")) +
  scale_fill_manual(values = c(baseline = "red3", vaccine = "bisque3"))

hosp_time_plot

ggsave(
  filename = glue("outputs/{sn}/figures/scenario_hospitalisations_over_time.png"),
  plot = hosp_time_plot,
  height = 5,
  width = 8
)

df <- inc_by_age_df |>
  group_by(prior_function_seed, age_group_plot, scenario) |>
  summarise(
    Incidence = sum(Incidence),
    Hospitalisations = sum(Hospitalisations),
    pop = mean(pop),
    .groups = "drop"
  ) |>
  group_by(age_group_plot, scenario) |>
  summarise(
    Incidence_lower = quantile(Incidence, probs = 0.025),
    Incidence_upper = quantile(Incidence, probs = 0.975),
    Incidence = median(Incidence),
    Hospitalisations_lower = quantile(Hospitalisations, probs = 0.025),
    Hospitalisations_upper = quantile(Hospitalisations, probs = 0.975),
    Hospitalisations = median(Hospitalisations),
    pop = median(pop),
    .groups = "drop"
  )

hosp_bar_plot <- ggplot(
  data = filter(df, age_group_plot %in% c("[60,65)", "[65,70)", "[70,75)", "[75, 90)")),
  aes(x = age_group_plot, y = Hospitalisations / pop * 100000, fill = scenario)
) +
  geom_col(position = "dodge") +
  geom_errorbar(
    aes(
      ymin = Hospitalisations_lower / pop * 100000,
      ymax = Hospitalisations_upper / pop * 100000
    ),
    width = 0.25,
    position = position_dodge(width = 0.9)
  ) +
  theme_minimal() +
  scale_fill_manual(values = c(baseline = "red3", vaccine = "bisque3")) +
  labs(x = "Age group", y = "Annual hospitalisations per 100,000")

hosp_bar_plot

ggsave(
  filename = glue("outputs/{sn}/figures/scenario_hospitalisations_by_age.png"),
  plot = hosp_bar_plot,
  height = 5,
  width = 8
)

inc_bar_plot <- ggplot(data = df, aes(x = age_group_plot, y = Incidence / pop * 1000, fill = scenario)) +
  geom_col(position = "dodge") +
  geom_errorbar(
    aes(
      ymin = Incidence_lower / pop * 1000,
      ymax = Incidence_upper / pop * 1000
    ),
    width = 0.25,
    position = position_dodge(width = 0.9)
  ) +
  theme_minimal() +
  scale_fill_manual(values = c(baseline = "red3", vaccine = "bisque3")) +
  labs(x = "Age group", y = "Annual infections per 1,000") +
  theme(axis.text.x = element_text(angle = 300, vjust = 0.5, hjust = 0))

inc_bar_plot

ggsave(
  filename = glue("outputs/{sn}/figures/scenario_infections_by_age.png"),
  plot = inc_bar_plot,
  height = 5,
  width = 10
)

hosp_age_time_plot <- ggplot(
  data = filter(inc_by_age_df_sum, age_group_plot %in% c("[60,65)", "[65,70)", "[70,75)", "[75, 90)")),
  aes(x = week, y = mH, colour = age_group_plot, fill = age_group_plot)
) +
  geom_line(linewidth = 0.9) +
  geom_ribbon(aes(ymin = lH, ymax = uH), alpha = 0.12, colour = NA) +
  facet_wrap(~scenario, ncol = 1) +
  theme_minimal() +
  labs(x = "Week", y = "Weekly hospitalisations", colour = "Age group", fill = "Age group")

hosp_age_time_plot

ggsave(
  filename = glue("outputs/{sn}/figures/scenario_hospitalisations_by_age_over_time.png"),
  plot = hosp_age_time_plot,
  height = 8,
  width = 9
)
