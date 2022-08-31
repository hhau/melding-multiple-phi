library(rjags)
library(coda)

source("scripts/common/mcmc-util.R")

model_data <- list(
  ti = 26,
  popcount = as.numeric(unlist(readRDS("rds/owls-example/count-data.rds")))
)

n_chain <- 6
n_iter <- 2e4

count_data_submodel <- jags.model(
  file = "scripts/owls-example/models/count-data.bug",
  data = model_data,
  n.chains = n_chain,
  n.adapt = n_iter / 2
)

parameters <- c(
  "phij",
  "phia",
  "fec",
  "im",  
  "v"
)

results <- coda.samples(
  model = count_data_submodel,
  variable.names = parameters,
  n.iter = n_iter
)

# reflow into iters X chains X parameters
results_array <- mcmc_list_to_array(
  results
)

saveRDS(
  object = results_array,
  file = "rds/owls-example/count-data-subposterior-samples.rds"
)
