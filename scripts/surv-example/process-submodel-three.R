library(tidyverse)
library(tidybayes)

source("scripts/common/mcmc-util.R")

submodel_three_output <- readRDS("rds/surv-example/submodel-three-output.rds")
samples <- submodel_three_output$samples

phi_23_names <- grep("^beta", names(samples[1, 1, ]), value = TRUE)
only_phi <- samples[, , phi_23_names]

point_est <- only_phi %>%
  array_to_mcmc_list() %>%
  gather_draws(beta[i, j]) %>%
  point_interval()

median_phi_stan_data <- list(
  long_beta_zero = point_est %>%
    filter(j == 1) %>%
    pull(.value),
  long_beta_one = point_est %>%
    filter(j == 2) %>%
    pull(.value)
)

saveRDS(
  file = "rds/surv-example/stage-one-phi-23-samples.rds",
  object = only_phi
)

saveRDS(
  file = "rds/surv-example/stage-one-phi-23-posterior-median.rds",
  object = median_phi_stan_data
)
