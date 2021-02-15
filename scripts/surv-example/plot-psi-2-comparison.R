library(tidyverse)
library(tidybayes)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

samples_melded <- readRDS("rds/surv-example/stage-two-psi-2-samples.rds")
samples_point <- readRDS("rds/surv-example/point-est-psi-2-samples.rds")

plot_tbl <- bind_rows(
  samples_melded %>%
    array_to_mcmc_list() %>%
    gather_draws(beta_zero, beta_one, hazard_gamma, alpha) %>%
    mutate(method = "melding"),
  samples_point %>%
    array_to_mcmc_list() %>%
    gather_draws(beta_zero, beta_one, hazard_gamma, alpha) %>%
    mutate(method = "point")
) %>%
  mutate(
    .variable = factor(
      x = .variable,
      levels = c("beta_zero", "beta_one", "hazard_gamma", "alpha"),
      labels = c("beta[0]", "beta[1]", "gamma", "alpha")
    ),
    method = as.factor(method)
  )

p1 <- ggplot(plot_tbl, aes(x = .value, lty = method)) +
  geom_density() +
  facet_wrap(vars(.variable), scales = "free", labeller = label_parsed) +
  labs(col = "Method") +
  scale_linetype_manual(
    values = c("melding" = "solid", "point" = "3313"),
    labels = c("melding" = "Melding", "point" = "Propagation")
  )

ggsave_halfheight(
  filename = "plots/surv-example/psi-2-method-comparison.pdf",
  plot = p1
)
