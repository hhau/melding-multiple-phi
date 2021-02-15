library(tidyverse)
library(tidybayes)

source("scripts/common/mcmc-util.R")

submodel_three_output <- readRDS("rds/surv-example/submodel-three-output.rds")
samples <- submodel_three_output$samples

phi_23_names <- grep("^beta_", names(samples[1, 1, ]), value = TRUE)
only_phi <- samples[, , phi_23_names]

point_est <- only_phi %>%
  array_to_mcmc_list() %>%
  gather_draws(beta_zero[i]) %>%
  point_interval()

median_phi_stan_data <- list(
  long_beta = point_est %>%
    filter(.variable == "beta_zero") %>%
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
