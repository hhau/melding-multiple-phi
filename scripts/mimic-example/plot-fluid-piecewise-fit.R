library(ggplot2)
library(tidybayes)
library(dplyr)

source('scripts/common/plot-settings.R')
source('scripts/common/setup-argparse.R')

parser$add_argument("--mimic-globals")
parser$add_argument("--cumulative-fluid-data")
parser$add_argument("--fluid-stan-data")
parser$add_argument("--fluid-plot-mu")
parser$add_argument("--output-plot-tbl")
args <- parser$parse_args()

source(args$mimic_globals)

cumulative_fluid_data <- readRDS(args$cumulative_fluid_data)
stan_data <- readRDS(args$fluid_stan_data)
plot_mu_samples <- readRDS(args$fluid_plot_mu)

plot_mu_tbl <- tibble(
  x = stan_data$x_mat_plot %>%
    t() %>%
    as.numeric(),
  i = rep(1 : stan_data$n_icu_stays, each = N_PLOT_POINTS),
  p = rep(1 : N_PLOT_POINTS, times = stan_data$n_icu_stays)
)

id_tbl <- tibble(
  i = 1 : stan_data$n_icu_stays,
  icustay_id = unique(cumulative_fluid_data$icustay_id)
)

plot_mu_interval <- plot_mu_samples %>%
  point_interval(.point = mean) %>%
  left_join(plot_mu_tbl, by = c('i', 'p')) %>%
  left_join(id_tbl, by = 'i')

base_plot <- ggplot(cumulative_fluid_data) +
  geom_point(aes(x = time_since_icu_adm, y = cumulative_value)) +
  facet_wrap(vars(icustay_id), scales = 'free')

p1 <- base_plot +
  geom_line(
    inherit.aes = FALSE,
    data = plot_mu_interval,
    mapping = aes(x = x, y = .value)
  ) +
  geom_ribbon(
    inherit.aes = FALSE,
    data = plot_mu_interval,
    mapping = aes(x = x, ymin = .lower, ymax = .upper),
    alpha = 0.3
  )

ggsave(
  filename = args$output,
  plot = p1,
  width = 18,
  height = 12
)

saveRDS(
  file = args$output_plot_tbl,
  object = plot_mu_interval
)

