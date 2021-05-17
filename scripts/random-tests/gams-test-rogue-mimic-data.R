# https://doi.org/10.1371/journal.pone.0226962 'data'
library(tidyverse)
library(mgcv)

train_data <- read_csv(
  file = "~/Downloads/pone.0226962.s002.csv"
)

max_stay_ids <- train_data %>% 
  group_by(icustay_id) %>%
  filter(first(pf) > 400) %>% 
  count(icustay_id, sort = TRUE) %>% 
  ungroup() %>% 
  slice((100 : 104) - 30)

plot_data <- train_data %>% 
  filter(icustay_id %in% max_stay_ids$icustay_id) %>% 
  select(X1, subject_id, pf, sf) %>% 
  pivot_longer(c(pf, sf), names_to = 'measurement_type', values_to = 'value') 

line_data <- plot_data %>% 
  group_by(subject_id, measurement_type) %>% 
  summarise(min_x = min(X1), max_x = max(X1)) %>% 
  filter(measurement_type == "pf")

ggplot(plot_data, aes(x = X1, y = value)) +
  geom_point() +
  facet_grid(cols = vars(subject_id), rows = vars(measurement_type), scales = 'free') +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    method.args = list(
      family = scat(min.df = 5)
    )
  ) +
  geom_hline(
    data = line_data,
    aes(yintercept = 300)
  )
