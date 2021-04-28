library(tidyverse)
library(tidybayes)
library(latex2exp)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

samples_melded <- readRDS("rds/mimic-example/stage-two-psi-2-samples.rds")
samples_point <- readRDS("rds/mimic-example/stage-two-median-inputs-psi-2-samples.rds")

n_theta <- samples_melded %>% 
  magrittr::extract(1, 1, ) %>% 
  names() %>% 
  grep(pattern = 'theta', x = .) %>% 
  length()

plot_tbl <- bind_rows(
  samples_melded %>%
    array_to_mcmc_list() %>%
    gather_draws(theta[e], hazard_gamma, alpha) %>%
    mutate(
      method = "melding",
      .value = ifelse(.variable == 'alpha', 1000 * .value, .value)
    ) %>% 
    unite('plot_var', c(.variable, e), na.rm = TRUE),
  samples_point %>%
    array_to_mcmc_list() %>%
    gather_draws(theta[e],  hazard_gamma, alpha) %>%
    mutate(
      method = "point",
      .value = ifelse(.variable == 'alpha', 1000 * .value, .value)
    ) %>% 
    unite('plot_var', c(.variable, e), na.rm = TRUE)
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
        "1000 %*% alpha"
      )
    ),
    method = as.factor(method)
  )

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

ggsave_fullpage(
  filename = "plots/mimic-example/psi-2-method-comparison.pdf",
  plot = p1
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
  aes(x = .value, colour = method)
) +
  geom_density() +
  facet_wrap(vars(.variable), scales = "free", labeller = label_parsed) +
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

ggsave_halfheight(
  filename = "plots/mimic-example/psi-2-method-comparison-small.pdf",
  plot = psmall
)

