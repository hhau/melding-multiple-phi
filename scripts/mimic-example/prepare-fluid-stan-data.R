library(dplyr)
library(purrr)
library(magrittr)

source("scripts/common/setup-argparse.R")
source("scripts/common/logger-setup.R")

parser$add_argument("--combined-pf-and-summarised-fluid-data")
parser$add_argument("--mimic-globals")
parser$add_argument("--output-cumulative-fluid")
args <- parser$parse_args()

source(args$mimic_globals)

summarised_data <- readRDS(args$combined_pf_and_summarised_fluid_data)

cumulative_fluid_data <- summarised_data %>%
  filter(value_type == 'fluids') %>%
  group_by(icustay_id) %>%
  arrange(icustay_id, time_since_icu_adm) %>%
  mutate(cumulative_value = cumsum(value))

saveRDS(
  file = args$output_cumulative_fluid,
  object = cumulative_fluid_data
)

n_icu_stays <- cumulative_fluid_data %>%
  pull(icustay_id) %>%
  n_distinct()

n_total_obs <- cumulative_fluid_data %>%
  nrow()

subset_vector <- cumulative_fluid_data %>%
  ungroup() %>%
  group_split(icustay_id) %>%
  map_int(~ nrow(.)) %>%
  cumsum() %>%
  add(1) %>%
  c(1, .)

breakpoint_limits <- cumulative_fluid_data %>%
  group_by(icustay_id) %>%
  summarise(
    lower = min(time_since_icu_adm),
    upper = max(time_since_icu_adm)
  )

emp_intercept_stats <- cumulative_fluid_data %>%
  group_by(icustay_id) %>%
  summarise(
    mean = mean(cumulative_value),
    sd = sd(cumulative_value)
  )

seq_vec <- Vectorize(seq.default, vectorize.args = c('from', 'to'))
x_mat_plot <- seq_vec(
  from = breakpoint_limits$lower,
  to = breakpoint_limits$upper,
  length.out = N_PLOT_POINTS
) %>% t()

stan_data <- list(
  n_icu_stays = n_icu_stays,
  n_total_obs = n_total_obs,
  n_plot_points = N_PLOT_POINTS,
  subset_vector = subset_vector,
  breakpoint_lower = breakpoint_limits$lower,
  breakpoint_upper = breakpoint_limits$upper,
  y_vec = cumulative_fluid_data$cumulative_value,
  x_vec = cumulative_fluid_data$time_since_icu_adm,
  x_mat_plot = x_mat_plot
)

saveRDS(
  file = args$output,
  object = stan_data
)
