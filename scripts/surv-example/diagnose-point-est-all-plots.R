library(tidyverse)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

point_est <- readRDS("rds/surv-example/point-est-psi-2-samples.rds")
point_est_1_meld_23 <- readRDS("rds/surv-example/point-est-1-meld-23-psi-2-samples.rds")
point_est_3_meld_12 <- readRDS("rds/surv-example/point-est-3-meld-12-psi-2-samples.rds")

par_names <- names(point_est[1, 1, ]) %>%
  grep("lp__", ., value = TRUE, invert = TRUE)

named_vec <- c(
  "beta_zero" = "theta[0]",
  'beta_one' = "theta[1]",
  "hazard_gamma" = "gamma",
  "alpha" = "alpha"
)

p_1 <- plot_worst_pars(point_est[, , par_names], named_vec)
p_2 <- plot_worst_pars(point_est_1_meld_23[, , par_names], named_vec)
p_3 <- plot_worst_pars(point_est_3_meld_12[, , par_names], named_vec)

ggsave_halfheight(
  filename = "plots/surv-example/point-est-diags.png",
  plot = p_1
)

ggsave_halfheight(
  filename = "plots/surv-example/point-est-1-meld-23-diags.png",
  plot = p_2
)

ggsave_halfheight(
  filename = "plots/surv-example/point-est-3-meld-12-diags.png",
  plot = p_3
)
