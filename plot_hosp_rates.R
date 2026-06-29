library(tidyverse)
df <- read_csv("data/hosp_rates2.csv") %>%
  mutate(age = factor(age, levels = c("0-6m", "6-12m", "12-24m", "24-60m", "5-9y","10-14y", "15-54y","55-59y","60-64y","65-69y","70-74y","75+y")
  ))

df <- df[1:12,]

ggplot(data = df, aes(x = age, y = rate)) +
  geom_bar(stat = "identity", fill = "darkorange3") +
  theme_minimal() +
  labs(x = "age", y = "rate per 100,000", title = "RSV hospitalisations") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 300, vjust = 0.5, hjust = 0))
