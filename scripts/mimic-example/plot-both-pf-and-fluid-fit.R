library(dplyr)
library(ggh4x)
library(tibble)
library(patchwork)

source('scripts/common/plot-settings.R')
source('scripts/common/setup-argparse.R')
source('scripts/mimic-example/GLOBALS.R')

parser$add_argument("--cumulative-fluid-data")
parser$add_argument("--fluid-plot-mu-tbl")
parser$add_argument("--pf-and-summarised-fluid-data")
parser$add_argument("--pf-plot-tbl")
parser$add_argument("--pf-event-time-samples-long")
parser$add_argument("--output-small")
args <- parser$parse_args()

fluid_plot_tbl <- readRDS(args$fluid_plot_mu_tbl) %>%
  select(-.variable) %>%
  mutate(value_type = 'fluids')

link_tbl <- fluid_plot_tbl %>%
  ungroup() %>%
  select(i, icustay_id) %>%
  distinct() %>%
  mutate(plot_id = i)

fluid_plot_tbl <- fluid_plot_tbl %>%
  left_join(link_tbl %>% select(-i), by = 'icustay_id')

cumulative_fluid_data <- readRDS(args$cumulative_fluid_data) %>%
  left_join(link_tbl, by = 'icustay_id')

pf_data <- readRDS(args$pf_and_summarised_fluid_data) %>%
  filter(value_type == 'pf') %>%
  left_join(link_tbl, by = 'icustay_id')

pf_plot_tbl <- readRDS(args$pf_plot_tbl) %>%
  left_join(link_tbl, by = 'i') %>%
  rename(.value = plot_mu) %>%
  mutate(value_type = 'pf')

event_time_tbl <- readRDS(args$pf_event_time_samples_long) %>%
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
    mapping = aes(x = time_since_icu_adm, y = cumulative_value),
    shape = 3,
    alpha = 0.75
  ) +
  geom_point(
    data = pf_data,
    mapping = aes(x = time_since_icu_adm, y = value)
  ) +
  geom_rug(
    inherit.aes = FALSE,
    data = event_time_tbl %>%
      filter(.variable == "event_time"),
    aes(x = .value),
    colour = highlight_col,
    alpha = 0.25
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
    ~ plot_id + value_type,
    scales = 'free',
    ncol = 4
  ) +
  xlab('Days since ICU admission') +
  theme(axis.title.y = element_blank())

ggsave(
  filename = args$output,
  plot = p1,
  width = 12,
  height = 48
)

p2 <- ggplot(
    interval_tbl %>%
      filter(plot_id %in% PLOT_IDS)
  ) +
  geom_point(
    data = pf_data %>%
      filter(plot_id %in% PLOT_IDS),
    mapping = aes(x = time_since_icu_adm, y = value)
  ) +
  geom_point(
    data = cumulative_fluid_data %>%
      filter(plot_id %in% PLOT_IDS),
    mapping = aes(x = time_since_icu_adm, y = cumulative_value),
    shape = 3,
    alpha = 0.75
  ) +
  geom_rug(
    data = event_time_tbl %>%
      filter(
        plot_id %in% PLOT_IDS,
        .variable == "event_time"
      ),
    aes(x = .value),
    colour = highlight_col,
    alpha = 0.01
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
    ~ plot_id + factor(
      value_type,
      levels = c('pf', 'fluids'),
      labels = c('P/F ratio (mmHg)', 'Cumulative fluid (L)')
    ),
    scales = 'free',
    ncol = 2
  ) +
  xlab(expression(italic(t) ~ (plain("Days since ICU admission")))) +
  theme(axis.title.y = element_blank())

ggsave_halfheight(
  filename = args$output_small,
  plot = p2,
)

patchwork_plot_list <- list()

for (ii in 1 : length(PLOT_IDS)) {
  is_right_column_plot <- ii == length(PLOT_IDS)

  if (is_right_column_plot) {
    pf_theme <- theme(axis.title = element_blank())
    fluid_theme <- theme(axis.title.y = element_blank())
  } else {
    pf_theme <- theme(axis.title.x = element_blank())
    fluid_theme <- NULL
  }

  # make the pf plot
  p_pf <- ggplot(
    interval_tbl %>%
      filter(
        plot_id %in% PLOT_IDS[ii],
        value_type == 'pf'
      )
  ) +
    geom_point(
      data = pf_data %>%
        filter(plot_id %in% PLOT_IDS[ii]),
      mapping = aes(x = time_since_icu_adm, y = value)
    ) +
    geom_rug(
      data = event_time_tbl %>%
        filter(
          plot_id %in% PLOT_IDS[ii],
          .variable == "event_time"
        ),
      aes(x = .value),
      colour = highlight_col,
      alpha = 0.01
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
    # xlab(expression(italic(t) ~ (plain("Days since ICU admission")))) +
    pf_theme +
    ylab('P/F ratio\n(mmHg)') +
    ggtitle(label = bquote(italic(i) == .(PLOT_IDS[ii])))

  # make the fluid plot
  p_fluid <- ggplot(
    interval_tbl %>%
      filter(
        plot_id %in% PLOT_IDS[ii],
        value_type == 'fluids'
      )
  ) +
    geom_point(
      data = cumulative_fluid_data %>%
        filter(plot_id %in% PLOT_IDS[ii]),
      mapping = aes(x = time_since_icu_adm, y = cumulative_value),
      shape = 3,
      alpha = 0.75
    ) +
    geom_line(
      mapping = aes(x = x, y = .value)
    ) +
    geom_ribbon(
      mapping = aes(x = x, ymin = .lower, ymax = .upper),
      alpha = 0.3
    ) +
    xlab(expression(italic(t) ~ (plain("Days since ICU admission")))) +
    ylab('Cumulative fluid\n(L)') +
    fluid_theme

  patchwork_plot_list[[((ii - 1) * 2) + 1]] <- p_pf
  patchwork_plot_list[[((ii - 1) * 2) + 2]] <- p_fluid
}

p3 <- wrap_plots(patchwork_plot_list, ncol = 2, nrow = 2, byrow = FALSE)
ggsave_halfheight(
  filename = 'plots/mimic-example/combined-pf-fluid-small-new.pdf',
  plot = p3,
  adjust_height = -2.5
)
