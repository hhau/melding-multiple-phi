library(dplyr)
library(ggh4x)
library(tibble)

source('scripts/common/plot-settings.R')

cumulative_fluid_data <- readRDS('rds/mimic-example/cumulative-summarised-fluid-data.rds')
fluid_plot_tbl <- readRDS('rds/mimic-example/fluid-data-piecewise-plot-mu-tbl.rds') %>% 
  select(-.variable) %>% 
  mutate(value_type = 'fluids')

link_tbl <- fluid_plot_tbl %>%
  ungroup() %>% 
  select(i, icustay_id) %>% 
  distinct()

pf_data <- readRDS('rds/mimic-example/combined-pf-and-summarised-fluids.rds') %>% 
  filter(value_type == 'pf')

pf_plot_tbl <- readRDS('rds/mimic-example/pf-data-and-bspline-plot-tbl.rds') %>% 
  left_join(link_tbl, by = 'i') %>% 
  rename(.value = plot_mu) %>% 
  mutate(value_type = 'pf')

event_time_tbl <- readRDS('rds/mimic-example/submodel-1-event-times-samples-long.rds') %>% 
  left_join(link_tbl, by = 'i') %>% 
  mutate(value_type = 'pf')

threshold_tbl <- tibble(
  yintercept = 300,
  value_type = 'pf'
)

interval_tbl <- bind_rows(fluid_plot_tbl, pf_plot_tbl)

p1 <- ggplot(interval_tbl) +
  geom_point(
    data = cumulative_fluid_data,
    mapping = aes(x = time_since_icu_adm, y = cumulative_value)
  ) +
  geom_point(
    data = pf_data,
    mapping = aes(x = time_since_icu_adm, y = value)
  ) +
  geom_rug(
    data = event_time_tbl,
    aes(x = .value),
    colour = highlight_col
  ) +
  geom_hline(
    data = threshold_tbl,
    mapping = aes(yintercept = yintercept)
  ) +
  geom_line(
    mapping = aes(x = x, y = .value)
  ) +
  geom_ribbon(
    mapping = aes(x = x, ymin = .lower, ymax = .upper),
    alpha = 0.3
  ) +
  facet_nested_wrap(
    ~ icustay_id + value_type,
    scales = 'free',
    ncol = 4
  )

ggsave(
  filename = 'plots/mimic-example/combined-pf-fluid-fit-plot.png',
  plot = p1,
  width = 12,
  height = 48
)

interesting_icustay_ids <- c(
  273811,
  268712,
  247952,
  213974
)

p2 <- ggplot(
    interval_tbl %>% 
      filter(icustay_id %in% interesting_icustay_ids)
  ) +
  geom_point(
    data = cumulative_fluid_data %>% 
      filter(icustay_id %in% interesting_icustay_ids),
    mapping = aes(x = time_since_icu_adm, y = cumulative_value)
  ) +
  geom_point(
    data = pf_data %>% 
      filter(icustay_id %in% interesting_icustay_ids),
    mapping = aes(x = time_since_icu_adm, y = value)
  ) +
  geom_rug(
    data = event_time_tbl %>% 
      filter(icustay_id %in% interesting_icustay_ids),
    aes(x = .value),
    colour = highlight_col
  ) +
  geom_hline(
    data = threshold_tbl,
    mapping = aes(yintercept = yintercept)
  ) +
  geom_line(
    mapping = aes(x = x, y = .value)
  ) +
  geom_ribbon(
    mapping = aes(x = x, ymin = .lower, ymax = .upper),
    alpha = 0.3
  ) +
  facet_nested_wrap(
    ~ icustay_id + value_type,
    scales = 'free',
    ncol = 4
  )

ggsave(
  filename = 'plots/mimic-example/combined-pf-fluid-fit-plot-small.png',
  plot = p2,
  width = 8,
  height = 6
)
