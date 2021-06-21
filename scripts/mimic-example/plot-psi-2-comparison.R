library(tidyverse)
library(tidybayes)
library(latex2exp)
library(sn)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")
source("scripts/common/setup-argparse.R")

parser$add_argument("--full-melding-poe-psi-2-samples")
parser$add_argument("--melding-poe-fix-phi-12-meld-phi-23-psi-2-samples")
parser$add_argument("--melding-poe-meld-phi-12-fix-phi-23-psi-2-samples")
parser$add_argument("--full-melding-logarthmic-psi-2-samples")
parser$add_argument("--melding-log-fix-phi-12-meld-phi-23-psi-2-samples")
parser$add_argument("--melding-log-meld-phi-12-fix-phi-23-psi-2-samples")
parser$add_argument("--both-fixed-psi-2-samples")
parser$add_argument("--median-event-time-data")
parser$add_argument("--output-small")
parser$add_argument("--output-alpha")
args <- parser$parse_args()

median_event_time_data <- readRDS(args$median_event_time_data)
samples_melded_poe <- readRDS(args$full_melding_poe_psi_2_samples)
samples_melded_poe_fix_phi_12_meld_phi_23 <- readRDS(args$melding_poe_fix_phi_12_meld_phi_23_psi_2_samples)
samples_melded_poe_meld_phi_12_fix_phi_23 <- readRDS(args$melding_poe_meld_phi_12_fix_phi_23_psi_2_samples)
samples_melded_log <- readRDS(args$full_melding_logarthmic_psi_2_samples)
samples_melded_log_fix_phi_12_meld_phi_23 <- readRDS(args$melding_log_fix_phi_12_meld_phi_23_psi_2_samples)
samples_melded_log_meld_phi_12_fix_phi_23 <- readRDS(args$melding_log_meld_phi_12_fix_phi_23_psi_2_samples)
samples_point <- readRDS(args$both_fixed_psi_2_samples)

n_samples <- dim(samples_melded_poe)[1]
n_chain <- dim(samples_melded_poe)[2]
n_to_draw <- n_samples * n_chain
n_theta <- samples_melded_poe %>%
  magrittr::extract(1, 1, ) %>%
  names() %>%
  grep(pattern = 'theta', x = .) %>%
  length()

log_crude_event_rate <- median_event_time_data[[2]] %>%
  pull()

samples_prior <- array(
  data = c(
    rnorm(n = n_to_draw, mean = log_crude_event_rate, sd = 0.5),
    rsn(n = (n_theta - 1) * n_to_draw, xi = 0, omega = 0.5, alpha = -1),
    rgamma(n_to_draw, shape = 9.05, rate = 8.72),
    rsn(n_to_draw, xi = 0, omega = 0.5, alpha = -2)
  ),
  dim = c(n_samples, n_chain, n_theta + 2),
  dimnames = dimnames(samples_melded_poe)
)

plot_tbl <- bind_rows(
  samples_prior %>%
    array_to_mcmc_list() %>%
    gather_draws(theta[b], hazard_gamma, alpha) %>%
    mutate(method = "prior", fix = "Not applic.") %>%
    unite('plot_var', c(.variable, b), na.rm = TRUE),
  samples_melded_poe %>%
    array_to_mcmc_list() %>%
    gather_draws(theta[b], hazard_gamma, alpha) %>%
    mutate(method = "melding-poe", fix = "none") %>%
    unite('plot_var', c(.variable, b), na.rm = TRUE),
  samples_melded_poe_fix_phi_12_meld_phi_23 %>%
    array_to_mcmc_list() %>%
    gather_draws(theta[b], hazard_gamma, alpha) %>%
    mutate(method = "melding-poe", fix = "phi_12") %>%
    unite('plot_var', c(.variable, b), na.rm = TRUE),
  samples_melded_poe_meld_phi_12_fix_phi_23 %>%
    array_to_mcmc_list() %>%
    gather_draws(theta[b], hazard_gamma, alpha) %>%
    mutate(method = "melding-poe", fix = "phi_23") %>%
    unite('plot_var', c(.variable, b), na.rm = TRUE),
  samples_melded_log %>%
    array_to_mcmc_list() %>%
    gather_draws(theta[b], hazard_gamma, alpha) %>%
    mutate(method = "melding-log", fix = "none") %>%
    unite('plot_var', c(.variable, b), na.rm = TRUE),
  samples_melded_log_fix_phi_12_meld_phi_23 %>%
    array_to_mcmc_list() %>%
    gather_draws(theta[b], hazard_gamma, alpha) %>%
    mutate(method = "melding-log", fix = "phi_12") %>%
    unite('plot_var', c(.variable, b), na.rm = TRUE),
  samples_melded_log_meld_phi_12_fix_phi_23 %>%
    array_to_mcmc_list() %>%
    gather_draws(theta[b], hazard_gamma, alpha) %>%
    mutate(method = "melding-log", fix = "phi_23") %>%
    unite('plot_var', c(.variable, b), na.rm = TRUE),
  samples_point %>%
    array_to_mcmc_list() %>%
    gather_draws(theta[b],  hazard_gamma, alpha) %>%
    mutate(method = "a_point", fix = "both") %>%
    unite('plot_var', c(.variable, b), na.rm = TRUE)
) %>%
  mutate(
    .variable = factor(
      x = plot_var,
      levels = c(
        sprintf('theta_%d', 1 : n_theta),
        "hazard_gamma",
        "alpha"
      ),
      labels = c(
        sprintf('theta[%d]', 1 : n_theta),
        "gamma",
        "alpha"
      )
    ),
    method = as.factor(method),
    fix = as.factor(fix)
  ) %>%
  filter(.iteration > 5)

p1 <- ggplot(
  plot_tbl,
  aes(x = .value, colour = method, linetype = fix)
) +
  geom_density() +
  facet_wrap(vars(.variable), scales = "free", labeller = label_parsed) +
  labs(colour = "Method", linetype = "Fixed") +
  scale_colour_manual(
    values = c(
      "prior" = 'grey',
      "melding-poe" = highlight_col,
      "melding-log" = greens[2],
      "a_point" = blues[1]
    ),
    labels = list(
      "prior" = TeX("Prior: $\\mathrm{p}_{2}(\\psi_{2})$"),
      "melding-poe" = "Chained melding, PoE pooling",
      "melding-log" = "Chained melding, log pooling",
      "a_point" = TeX("Fix $\\phi_{1 \\bigcap 2}$ and $\\phi_{2 \\bigcap 3}$")
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Not applic." = "solid",
      "none" = "solid",
      "both" = "solid",
      "phi_12" = "dashed",
      "phi_23" = "dotted"
    ),
    labels = list(
      "Not applic." = "Not applic.",
      "none" = "Full chained melding",
      "both" = TeX("Fix $\\phi_{1 \\bigcap 2}$ and $\\phi_{2 \\bigcap 3}$"),
      "phi_12" = TeX("Fix $\\phi_{1 \\bigcap 2}$, meld $\\phi_{2 \\bigcap 3}$"),
      "phi_23" = TeX("Meld $\\phi_{1 \\bigcap 2}$, fix $\\phi_{2 \\bigcap 3}$")
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
  'theta_3',
  'theta_17'
)

psmall <- ggplot(
  plot_tbl %>%
    filter(plot_var %in% interesting_subplots),
  aes(x = .value, colour = method, linetype = fix)
) +
  geom_density() +
  facet_wrap(
    vars(.variable),
    scales = "free",
    labeller = label_parsed,
    ncol = 2
  ) +
  labs(colour = "Method", linetype = "Fixed") +
  scale_colour_manual(
    values = c(
      "prior" = 'grey',
      "melding-poe" = highlight_col,
      "melding-log" = greens[2],
      "a_point" = blues[1]
    ),
    labels = list(
      "prior" = TeX("Prior: $\\mathrm{p}_{2}(\\psi_{2})$"),
      "melding-poe" = "Chained melding, PoE pooling",
      "melding-log" = "Chained melding, log pooling",
      "a_point" = TeX("Fix $\\phi_{1 \\bigcap 2}$ and $\\phi_{2 \\bigcap 3}$")
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Not applic." = "solid",
      "none" = "solid",
      "both" = "solid",
      "phi_12" = "dashed",
      "phi_23" = "dotted"
    ),
    labels = list(
      "Not applic." = "Not applic.",
      "none" = "Full chained melding",
      "both" = TeX("Fix $\\phi_{1 \\bigcap 2}$ and $\\phi_{2 \\bigcap 3}$"),
      "phi_12" = TeX("Fix $\\phi_{1 \\bigcap 2}$, meld $\\phi_{2 \\bigcap 3}$"),
      "phi_23" = TeX("Meld $\\phi_{1 \\bigcap 2}$, fix $\\phi_{2 \\bigcap 3}$")
    )
  ) +
  theme(axis.title = element_blank())

ggsave_halfheight(
  filename = args$output_small,
  plot = psmall
)

interval_tbl <- plot_tbl %>%
  filter(plot_var == 'alpha') %>%
  group_by(method, fix, plot_var) %>%
  tidybayes::point_interval(
    .value,
    .width = c(0.5, 0.8, 0.95, 0.99)
  ) %>%
  mutate(plot_y = interaction(method, fix))

label_tbl <- tribble(
  ~plot_y, ~lab_y, ~val_y,
  "prior.Not applic.", "a_prior", 11,
  "melding-poe.none", "b_poe_none", 9,
  "melding-poe.phi_12", "bb_poe_phi_12", 8,
  "melding-poe.phi_23", "bbb_poe_phi_23", 7,
  "melding-log.none", "c_log_none", 5,
  "melding-log.phi_12", "cc_log_phi_12", 4,
  "melding-log.phi_23", "ccc_log_phi_23", 3,
  "a_point.both", "d_point_both", 1
)

plot_two_tbl <- interval_tbl %>%
  left_join(label_tbl, by = 'plot_y') %>%
  arrange(desc(val_y)) %>%
  ungroup() %>%
  mutate(
    val_y = factor(
      x = val_y,
      levels = c(11, 9, 8, 7, 5, 4, 3, 1),
      labels = c(
          TeX("Prior: NA"),
          TeX("PoE: none"),
          TeX("PoE: fix $\\phi_{1 \\bigcap 2}$"),
          TeX("PoE: fix $\\phi_{2 \\bigcap 3}$"),
          TeX("Log: none"),
          TeX("Log: fix $\\phi_{1 \\bigcap 2}$"),
          TeX("Log: fix $\\phi_{2 \\bigcap 3}$"),
          TeX("Point: both")
        )
    )
  )

p_alpha <- ggplot(
  data = plot_two_tbl,
  aes(
    x = .value,
    xmin = .lower,
    xmax = .upper,
    y = val_y,
    col = method,
    alpha = 1 - .width
  )
) +
  geom_interval(size = rel(5)) +
  geom_point(
    aes(x = .value, y = val_y),
    pch = '|',
    size = rel(6),
    show.legend = FALSE,
    stroke = 2
  ) +
  scale_alpha(guide = FALSE) +
  theme(
    axis.title.y = element_blank()
  ) +
  xlab(expression(alpha)) +
  scale_y_discrete(
    limits = rev,
    labels = function(x) parse(text = x)
  ) +
  scale_colour_manual(
    values = c(
      "prior" = 'grey',
      "melding-poe" = highlight_col,
      "melding-log" = greens[2],
      "a_point" = 'black'
    ),
    labels = list(
      "prior" = TeX("Prior: $\\mathrm{p}_{2}(\\psi_{2})$"),
      "melding-poe" = "PoE pooling",
      "melding-log" = "Log pooling",
      "a_point" = TeX("Fix $\\phi_{1 \\bigcap 2}$ and $\\phi_{2 \\bigcap 3}$")
    ),
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(colour = "")

ggsave_halfheight(
  plot = p_alpha,
  filename = args$output_alpha
)
