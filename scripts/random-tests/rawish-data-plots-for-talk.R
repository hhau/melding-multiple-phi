library(tidyverse)

source('scripts/common/plot-settings.R')

both_data <- readRDS('rds/mimic-example/combined-pf-and-summarised-fluids.rds')
pf_data <- both_data %>% 
  filter(value_type == 'pf')

fluid_data <- readRDS('rds/mimic-example/cumulative-summarised-fluid-data.rds') 
  
plot_icustay_ids <- c(
  202327,
  252051,
  272300
)

p1_small <- ggplot(
  pf_data %>% filter(icustay_id %in% plot_icustay_ids), 
  aes(x = time_since_icu_adm, y = value)
) + 
  geom_point() +
  geom_hline(yintercept = 300, lty = 'dashed', col = highlight_col) +
  facet_grid(cols = vars(icustay_id), scales = 'free_x') +
  xlab('Days since ICU admission') +
  ylab('P/F ratio')

ggsave(
  plot = p1_small,
  filename = 'plots/mimic-example/pf-data-raw-small.pdf',
  width = 4.25 * 1.33,
  height = 2 * 1.33
)

p2_small <- ggplot(
  fluid_data %>% filter(icustay_id %in% plot_icustay_ids), 
  aes(x = time_since_icu_adm, y = cumulative_value)
) + 
  geom_point() +
  facet_grid(cols = vars(icustay_id), scales = 'free_x') +
  xlab('Days since ICU admission') +
  ylab('Cumulative fluid balance')

ggsave(
  plot = p2_small,
  filename = 'plots/mimic-example/fluid-data-cumulative-small.pdf',
  width = 4.25 * 1.33,
  height = 2 * 1.33
)
