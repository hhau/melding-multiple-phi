library(tidyverse)
library(ggh4x)
library(mgcv)

pf_data <- readRDS('rds/mimic-tests/pf-cohort-and-measurements.rds')
fluid_data <- readRDS('rds/mimic-tests/cumulative-fluids-all-patients.rds')

pf_time_extrema <- pf_data %>% 
  group_by(icustay_id) %>% 
  summarise(
    pf_min_time = min(time_since_icu_adm),
    pf_max_time = max(time_since_icu_adm)
  ) 

fluid_time_extrema <- fluid_data %>%
  filter(icustay_id %in% unique(pf_data$icustay_id)) %>% 
  group_by(icustay_id) %>% 
  summarise(
    fluid_min_time = min(time_since_icu_adm),
    fluid_max_time = max(time_since_icu_adm)
  )

all_data <- inner_join(
  pf_time_extrema,
  fluid_time_extrema,
  by = 'icustay_id'
) %>% 
  arrange(icustay_id) %>% 
  filter(complete.cases(.)) %>% 
  rowwise() %>% 
  mutate(
    overlap = max(
      0,
      min(pf_max_time, fluid_max_time) - max(pf_min_time, fluid_min_time)
    ),
    overlap_as_proportion_of_possible = overlap / (
      max(pf_max_time, fluid_max_time) - min(pf_min_time, fluid_min_time)
    )
  )

overlap_cohort <- all_data %>% 
  filter(overlap_as_proportion_of_possible > 0.9)

plot_tbl <- pf_data %>% 
  filter(icustay_id %in% overlap_cohort$icustay_id) %>% 
  select(
    -c(subject_id, hadm_id, charttime, specimen_pred, intime),
    value = pf
  ) %>% 
  mutate(value_type = 'pf')

both_plot_tbl <- fluid_data %>% 
  filter(icustay_id %in% overlap_cohort$icustay_id) %>% 
  select(-c(origin, hadm_id, cumulative_fluids)) %>% 
  mutate(value_type = 'fluids') %>% 
  bind_rows(plot_tbl)

p2 <- ggplot(both_plot_tbl, aes(x = time_since_icu_adm, y = value)) +
  geom_point() +
  facet_nested_wrap(vars(icustay_id, value_type), scales = 'free_y') +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    method.args = list(
      family = scat(min.df = 5)
    )
  ) +
  geom_hline(
    data = tibble(
      value_type = 'pf',
      yintercept = 300
    ),
    mapping = aes(
      yintercept = yintercept
    )
  ) +
  xlab('Days since start of this ICU admission') +
  theme(
    strip.text = element_text(margin = margin(t = 2, b = 2, l = 50, r = 50))
  )

ggsave(
  filename = "plots/mimic-tests/reasonable-pf-and-fluids-data.pdf",
  plot = p2,
  height = 30,
  width = 45,
  units = 'in'
)
