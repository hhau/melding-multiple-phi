library(tidyverse)
library(tidybayes)

source("scripts/common/mcmc-util.R")

submodel_one_output <- readRDS("rds/surv-example/submodel-one-output.rds")
samples <- submodel_one_output$samples

phi_12_names <- grep("^event", names(samples[1, 1, ]), value = TRUE)
only_phi <- samples[, , phi_12_names]

point_est <- only_phi %>%
  array_to_mcmc_list() %>%
  gather_draws(event_time[i], event_indicator[i]) %>%
  point_interval()

median_phi_stan_data <- list(
  event_indicator = point_est %>%
    filter(.variable == "event_indicator") %>%
    pull(.value),
  event_time = point_est %>%
    filter(.variable == "event_time") %>%
    pull(.value)
)

saveRDS(
  file = "rds/surv-example/stage-one-phi-12-samples.rds",
  object = only_phi
)

saveRDS(
  file = "rds/surv-example/stage-one-phi-12-posterior-median.rds",
  object = median_phi_stan_data
)
