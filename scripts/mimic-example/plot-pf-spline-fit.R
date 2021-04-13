library(tidybayes)
library(ggdist)
library(dplyr)
library(abind)

source('scripts/common/plot-settings.R')
source("scripts/common/setup-argparse.R")
source("scripts/common/logger-setup.R")

bind_named_sublists <- function(outer_list, name, along_dim) {
  lapply(outer_list, function(sub_list) sub_list[[name]]) %>%
    abind(along = along_dim)
}

parser$add_argument("--mimic-pf-plot-mu")
parser$add_argument("--mimic-pf-data-list")
parser$add_argument("--mimic-pf-event-time-long")
parser$add_argument("--mimic-globals")
parser$add_argument("--combined-pf-and-summarised-fluid-data")
args <- parser$parse_args()

source(args$mimic_globals)

list_data <- readRDS(args$mimic_pf_data_list)
post_samples <- readRDS(args$mimic_pf_plot_mu)
event_time_samples <- readRDS(args$mimic_pf_event_time_long)

n_icu_stays <- length(list_data)
pf_data <- readRDS(args$combined_pf_and_summarised_fluid_data) %>%
  filter(value_type == 'pf')

plot_tbl <- post_samples %>%
  mean_qi()

plot_x_df <- bind_named_sublists(list_data, 'x_plot_seq', 1) %>%
  as_tibble() %>%
  mutate(
    i = rep(1 : n_icu_stays, each = N_PLOT_POINTS),
    p = rep(1 : N_PLOT_POINTS, times = n_icu_stays)
  )

plot_tbl <- plot_tbl %>%
  left_join(plot_x_df) %>%
  rename(x = value)

point_tbl <- pf_data %>%
  mutate(i = as.numeric(as.factor(icustay_id)))

p1 <- ggplot(data = plot_tbl) +
  geom_line(aes(x = x, y = plot_mu)) +
  geom_ribbon(aes(x = x, ymin = .lower, ymax = .upper), alpha = 0.2) +
  geom_point(
    data = point_tbl,
    inherit.aes = FALSE,
    mapping = aes(x = time_since_icu_adm, y = value)
  ) +
  geom_hline(
    data = tibble(
      yintercept = 300
    ),
    mapping = aes(yintercept = yintercept)
  ) +
  geom_rug(
    data = event_time_samples %>%
      filter(.variable == "event_time"),
    inherit.aes = FALSE,
    mapping = aes(x = .value),
    alpha = 0.5
  ) +
  facet_wrap(vars(i), scales = 'free')

ggsave(
  filename = args$output,
  plot = p1,
  width = 30,
  height = 20
)