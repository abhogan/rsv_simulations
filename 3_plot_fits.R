rm(list = ls())
source("data/data.r")
source("R/set_up_variables.r")

sn <- "alpha_fit_nu200_omega1"
# contact matrix plot
nAges <- length(contact_population_list$age_limits)

matrix_plot <- contact_population_list$matrix_per_person |> as.data.frame() |>
  mutate(age_group = rownames(contact_population_list$matrix_per_person)) |>
         pivot_longer(cols = c(-age_group), values_to = "contacts", names_to = "contact_age_group")

matrix_plot$age_group <- factor(matrix_plot$age_group, levels = rownames(contact_population_list$matrix_per_person))
matrix_plot$contact_age_group <- factor(matrix_plot$contact_age_group, levels = rownames(contact_population_list$matrix_per_person))

ggsave(
  ggplot(data = matrix_plot,
       aes(x = age_group, y = contact_age_group, fill = contacts)) +
         geom_tile() +
    scale_fill_gradientn(colors = rev(heat.colors(70))) +
    ylab("Contact age group (years)") + xlab("Age group (years)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 15)),
  file = "outputs/new/figures/contact_matrix_heatmap.pdf",
  device = pdf,
  height = 10,
  width = 11
  )

###############################################
##### reading in the runs with model fits #####
###############################################

fit_all <- readRDS(file = glue("outputs/{sn}/model_fit_ABC_SMC_100_min.rds"))
G <- length(fit_all$fitted_parameters)

fit <- fit_all$fitted_parameters[[G]] |>
  as.data.frame()
colnames(fit) <- c("b0", "b1", "phi", "alpha")

model_runs <- readRDS(glue("outputs/{sn}/model_runs_100_min.rds"))
inc_df <- readRDS(glue("outputs/{sn}/inc_df.rds"))
inc_by_age_df <- readRDS(glue("outputs/{sn}/inc_by_age_df.rds"))
inc_by_age_df_sum <- readRDS(glue("outputs/{sn}/inc_by_age_df_sum.rds"))
max_week_df <- readRDS(glue("outputs/{sn}/max_week_df.rds"))
ar_df <- readRDS(glue("outputs/{sn}/ar_df.rds"))
#ar_df_c <- readRDS(glue("outputs/{sn}/ar_df_c.rds"))
ar_df_03 <- readRDS(glue("outputs/{sn}/ar_df_03.rds")) %>%
  mutate(age = 1:3) %>%
  mutate(age = factor(age, labels = c("0-1y", "0-2y", "0-3y")))

ar_df_65 <- ar_df[14:16,] %>%
  summarise(ar_median = mean(mAR),
            ar_lower = mean(lAR),
            ar_upper = mean(uAR)) %>%
  mutate(age = "65+y")

ar_df_03 <- rbind(ar_df_03, ar_df_65)

#####################
##### incidence #####
#####################

age_inc_plot <- ggplot(data = inc_by_age_df_sum |> mutate(age_group_plot = factor(age_group_plot)),
       aes(x = week, y = mI, ymin = lI, ymax = uI, group = factor(age_group_plot))) +
  geom_ribbon(alpha = 0.1, aes(fill = age_group_plot)) +
  geom_line(aes(col = age_group_plot), linewidth = 1) +
  theme_bw()  + theme(text = element_text(size = 17)) +
  scale_colour_viridis_d(name = "Age group") +
  scale_fill_viridis_d(name = "Age group") +
  ylab("Simulated weekly incidence per 1000 people") +
  xlab("Week of year")

age_inc_plot

peak_inc_plot <- ggplot(data = max_week_df, aes(x = age_group_plot, ymin = lW, y = mW, ymax = uW)) +
  geom_pointrange(size = 0.75) +
  scale_fill_viridis_d(name = "Age group (years)") +
  theme_bw() + theme(text = element_text(size = 17), legend.position = "none") +
  xlab("Age group") + ylab("Week of peak incidence") +
  #scale_y_continuous(limits = c(0, 50)) +
  theme(axis.text.x = element_text(angle = 300, vjust = 0.5, hjust=0))

peak_inc_plot

inc_plot <- ggplot(data = inc_df, aes(x = week, y = median_sc, ymin = lower_sc, ymax = upper_sc)) +
  geom_bar(inherit.aes = FALSE,
           data = data.frame(time = case_times_pia,
                             inc = cases_pia / max(cases_pia)),
           aes(x = time, y = inc), stat = "identity", alpha = 0.35,
           col = "grey50") +
  geom_ribbon(alpha = 0.5, fill = "skyblue") +
  geom_line(linewidth = 1, col = "skyblue") +
  #geom_point(col = "skyblue") +
  theme_bw() + theme(text = element_text(size = 17)) +
  xlab("Weeks") + ylab("Scaled incidence")

inc_plot

ggsave(glue("outputs/{sn}/figures/inc_plot.png"), height = 5, width = 8)

inc_plot <- ggplot(data = inc_df, aes(x = week, y = median_tot/27536874*100000, ymin = lower_tot/27536874*100000, ymax = upper_tot/27536874*100000)) +
  geom_ribbon(alpha = 0.5, fill = "skyblue") +
  geom_line(linewidth = 1, col = "skyblue") +
  #geom_point(col = "skyblue") +
  theme_bw() + theme(text = element_text(size = 17)) +
  xlab("week") + ylab("incidence per 100,000")

inc_plot

ggsave(glue("outputs/{sn}/figures/inc_plot_all.png"), height = 5, width = 8)


########################
##### attack rates #####
########################

annual_ar_plot <- ggplot(data = ar_df, aes(x = age_group_plot, y = mAR, ymin = lAR, ymax = uAR)) +
  geom_pointrange(size = 0.75) +
  ylab("modelled annual attack rate") +
  xlab("age group (years)") +
  theme_bw() + theme(text = element_text(size = 17), legend.position = "none") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  theme(axis.text.x = element_text(angle = 300, vjust = 0.5, hjust=0))

annual_ar_plot

ggsave(glue("outputs/{sn}/figures/annual_ar_plot.png"), height = 5, width = 8)

ar_plot <- ggplot(data = ar_df_c,
       aes(x = age_continuous, y = ar_m, ymin = ar_l, ymax = ar_u,
           fill = factor(cohort_birth), group = cohort_birth)) +
  theme_bw() + theme(text = element_text(size = 17)) +
  #geom_ribbon(alpha = 0.1) +
  geom_line(aes(col = factor(cohort_birth)), linewidth = 0.5) +
  ylab("Cumulative attack rate (primary infection)") +
  coord_cartesian(xlim = c(0, 3)) +
  geom_segment(data = data.frame(xmin = ages_ar - 1, xmax = ages_ar, ar = attack_rate, study = ref_ar),
               inherit.aes = FALSE,
               aes(x = xmin, y = ar, yend = ar, xend = xmax, linetype = study)) +
  guides(col = "none", fill = "none") +
  xlab("Age (years)") + scale_y_continuous(labels = scales::percent)#labs(tag = "C")

ar_plot

ages_ar_plot <- factor(ages_ar, levels = c(1,2,3,70), labels = c("0-1y", "0-2y", "0-3y", "65+y"))
ar_data <- data.frame(age = ages_ar_plot, ar = attack_rate, study = ref_ar, ar_lower = NA, ar_upper = NA)
ar_03_plot <- ggplot(data = ar_df_03,
                     aes(x = age, y = ar_median, ymin = ar_lower, ymax = ar_upper)) +
  geom_pointrange() +
  theme_bw() + theme(text = element_text(size = 17)) +
  ylab("attack rate\n(primary infection)") +
  geom_jitter(data = ar_data,
             aes(x = age, y = ar, col = study, fill = study, shape = study), size = 3) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) #+
  #theme(axis.text.x = element_text(angle = 300, vjust = 0.5, hjust=0))

ar_03_plot
ggsave(glue("outputs/{sn}/figures/annual_ar_younger.png"), height = 4, width = 6)

########################
# 
# ggsave(file = "figures/calibration_plot.pdf",
#        plot = (ar_plot / inc_plot) + plot_annotation(tag_levels = "A"),
#        units = "cm",
#        height = 30,
#        width = 37.5,
#        dpi = 300
#        )

#################
##### plots #####
#################

# model calibration

posterior_plot <- (ggplot(data = fit,
                         aes(x = b0)) +
  geom_density(linewidth = 1, fill = "skyblue", alpha = 0.5) + theme_bw() +
  theme(text = element_text(size = 17)) +
  scale_y_continuous(breaks = seq(0, 250, 50)) +
  xlab(bquote(Transmission~rate~coefficient~(b[0]))) +
  ylab("Posterior density") + labs(tag = "A") +
  ggplot(data = fit,
         aes(x = b1)) +
  xlab(bquote(Seasonal~forcing~amplitude~(b[1]))) +
  ylab("Posterior density") +
  geom_density(linewidth = 1, fill = "skyblue", alpha = 0.5) + theme_bw() +
  theme(text = element_text(size = 17)) +
  ggplot(data = fit,
         aes(x = phi)) +
  geom_density(linewidth = 1, fill = "skyblue", alpha = 0.5) + theme_bw() +
  theme(text = element_text(size = 17)) +
  scale_x_continuous(breaks = seq(125, 365, 50)) +
  xlab(expression(paste("Seasonal forcing phase shift (", phi, ")"))) +
  ylab("Posterior density") +
  ggplot(data = fit,
         aes(x = alpha)) +
    geom_density(linewidth = 1, fill = "skyblue", alpha = 0.5) + theme_bw() +
    theme(text = element_text(size = 17)) +
    scale_x_continuous(breaks = seq(alpha_lower, alpha_upper, 0.1)) +
    xlab(expression(paste("Age-independent susceptibility factor (", alpha, ")"))) +
    ylab("Posterior density") +
  plot_layout(axis_titles = "collect"))
posterior_plot

labels <- c(expression(paste("Transmission rate coefficient (",b[0],")")),
            expression(paste("Seasonal forcing amplitude (",b[1],")")),
            expression(paste("Phase-shift (",phi,")")),
            expression(paste("Age-independent susceptibility factor (",alpha,")"))
            )

labels_fixed <- c("Transmission~rate~coefficient~(b[0])",
            "Seasonal~forcing~amplitude~(b[1])",
            "Phase-shift~(phi)",
            "Age-independent~susceptibility~factor~(alpha)")

smooth_fn <- function(data, mapping, method = "lm", ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.5, col = "grey50") +
    geom_smooth(method = method, ...)
}

pairs_plot <- ggpairs(data = fit[,c("b0", "b1", "phi", "alpha")],
        lower = list(continuous = wrap(smooth_fn, method = "loess", color = "black", se = FALSE)),
        diag = list(continuous = wrap("densityDiag", alpha = 0.5, fill = "skyblue")),
        upper =  NULL,
        #facet = list(labeller = custom_labeler),
        columnLabels = labels_fixed,
        labeller = "label_parsed"
        ) +
  theme_bw() +
  theme(text = element_text(size = 17),
        strip.text = element_text(size = 10))

ggsave(
  filename = "figures/model_calibration.pdf",
  plot = pairs_plot,
  units = "cm",
  height = 40,
  width = 40,
  dpi = 300
)


