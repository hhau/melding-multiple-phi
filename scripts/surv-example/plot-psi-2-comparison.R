library(tidyverse)
library(tidybayes)
library(latex2exp)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

samples_melded <- readRDS("rds/surv-example/stage-two-psi-2-samples.rds")
samples_point <- readRDS("rds/surv-example/point-est-psi-2-samples.rds")
samples_point_1_melded_23 <- readRDS("rds/surv-example/point-est-1-meld-23-psi-2-samples.rds")
samples_point_3_melded_12 <- readRDS("rds/surv-example/point-est-3-meld-12-psi-2-samples.rds")

plot_tbl <- bind_rows(
  samples_melded %>%
    array_to_mcmc_list() %>%
    gather_draws(theta_zero, theta_one, hazard_gamma, alpha) %>%
    mutate(method = "melding"),
  samples_point %>%
    array_to_mcmc_list() %>%
    gather_draws(theta_zero, theta_one, hazard_gamma, alpha) %>%
    mutate(method = "point"),
  samples_point_1_melded_23 %>%
    array_to_mcmc_list() %>%
    gather_draws(theta_zero, theta_one, hazard_gamma, alpha) %>%
    mutate(method = "point-1-meld-23"),
  samples_point_3_melded_12 %>%
    array_to_mcmc_list() %>%
    gather_draws(theta_zero, theta_one, hazard_gamma, alpha) %>%
    mutate(method = "point-3-meld-12")
) %>%
  mutate(
    .variable = factor(
      x = .variable,
      levels = c("theta_zero", "theta_one", "hazard_gamma", "alpha"),
      labels = c("theta[0]", "theta[1]", "gamma", "alpha")
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
      "point" = blues[1],
      "point-1-meld-23" = blues[2],
      "point-3-meld-12" = blues[3]
    ),
    labels = list(
      "melding" = "Chained melding",
      "point" = TeX("Fix $\\phi_{1 \\bigcap 2}$ and $\\phi_{2 \\bigcap 3}$"),
      "point-1-meld-23" = TeX("Fix $\\phi_{1 \\bigcap 2}$, meld $\\phi_{2 \\bigcap 3}$"),
      "point-3-meld-12" = TeX("Meld $\\phi_{1 \\bigcap 2}$, fix $\\phi_{2 \\bigcap 3}$")
    )
  ) +
  scale_linetype_manual(
    values = c(
      "melding" = "solid",
      "point" = "solid",
      "point-1-meld-23" = "3111",
      "point-3-meld-12" = "3111"
    ),
    labels = list(
      "melding" = "Chained melding",
      "point" = TeX("Fix $\\phi_{1 \\bigcap 2}$ and $\\phi_{2 \\bigcap 3}$"),
      "point-1-meld-23" = TeX("Fix $\\phi_{1 \\bigcap 2}$, meld $\\phi_{2 \\bigcap 3}$"),
      "point-3-meld-12" = TeX("Meld $\\phi_{1 \\bigcap 2}$, fix $\\phi_{2 \\bigcap 3}$")
    )
  ) +
  theme(axis.title = element_blank())

ggsave_halfheight(
  filename = "plots/surv-example/psi-2-method-comparison.pdf",
  plot = p1
)
