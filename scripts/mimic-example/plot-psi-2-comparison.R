library(tidyverse)
library(tidybayes)
library(latex2exp)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")
source("scripts/common/setup-argparse.R")

parser$add_argument("--full-melding-psi-2-samples")
parser$add_argument("--both-fixed-psi-2-samples")
parser$add_argument("--output-small")
args <- parser$parse_args()

samples_melded <- readRDS(args$full_melding_psi_2_samples)
samples_point <- readRDS(args$both_fixed_psi_2_samples)

n_theta <- samples_melded %>%
  magrittr::extract(1, 1, ) %>%
  names() %>%
  grep(pattern = 'theta', x = .) %>%
  length()

plot_tbl <- bind_rows(
  samples_melded %>%
    array_to_mcmc_list() %>%
    gather_draws(theta[b], hazard_gamma, dd_gamma, alpha) %>%
    mutate(
      method = "melding",
      .value = ifelse(.variable == 'alpha', 1000 * .value, .value)
    ) %>%
    unite('plot_var', c(.variable, b), na.rm = TRUE),
  samples_point %>%
    array_to_mcmc_list() %>%
    gather_draws(theta[b],  hazard_gamma, dd_gamma, alpha) %>%
    mutate(
      method = "point",
      .value = ifelse(.variable == 'alpha', 1000 * .value, .value)
    ) %>%
    unite('plot_var', c(.variable, b), na.rm = TRUE)
) %>%
  mutate(
    .variable = factor(
      x = plot_var,
      levels = c(
        sprintf('theta_%d', 1 : n_theta),
        "hazard_gamma",
        "dd_gamma",
        "alpha"
      ),
      labels = c(
        sprintf('theta[%d]', 1 : n_theta),
        "gamma[1]",
        "gamma[2]",
        "1000 %*% alpha"
      )
    ),
    method = as.factor(method)
  ) %>%
  filter(.iteration > 5)

p1 <- ggplot(
  plot_tbl,
  aes(x = .value, colour = method, linetype = method)
) +
  geom_density() +
  facet_wrap(vars(.variable), scales = "free", labeller = label_parsed) +
  labs(colour = "Method", linetype = "Method") +
  scale_colour_manual(
    values = c(
      "melding" = highlight_col,
      "point" = blues[1]
    ),
    labels = list(
      "melding" = "Chained melding",
      "point" = TeX("Fix $\\phi_{1 \\bigcap 2}$ and $\\phi_{2 \\bigcap 3}$")
    )
  ) +
  scale_linetype_manual(
    values = c(
      "melding" = "solid",
      "point" = "solid"
    ),
    labels = list(
      "melding" = "Chained melding",
      "point" = TeX("Fix $\\phi_{1 \\bigcap 2}$ and $\\phi_{2 \\bigcap 3}$")
    )
  ) +
  theme(axis.title = element_blank())

ggsave_base(
  filename = args$output,
  plot = p1,
  height = display_settings$full_page_plot_height,
  width = display_settings$full_page_plot_width + 10
)

interesting_subplots <- c(
  'alpha',
  'hazard_gamma',
  'dd_gamma',
  'theta_3',
  'theta_17'
)

psmall <- ggplot(
  plot_tbl %>%
    filter(plot_var %in% interesting_subplots),
  aes(x = .value, colour = method)
) +
  geom_density() +
  facet_wrap(
    vars(.variable),
    scales = "free",
    labeller = label_parsed,
    ncol = 2
  ) +
  labs(colour = "Method") +
  scale_colour_manual(
    values = c(
      "melding" = highlight_col,
      "point" = blues[1]
    ),
    labels = list(
      "melding" = "Chained melding",
      "point" = TeX("Fix $\\phi_{1 \\bigcap 2}$ and $\\phi_{2 \\bigcap 3}$")
    )
  ) +
  theme(axis.title = element_blank())

ggsave_fullpage(
  filename = args$output_small,
  plot = psmall
)
