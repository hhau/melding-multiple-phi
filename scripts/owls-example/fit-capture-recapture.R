library(rjags)
library(coda)

source("scripts/common/mcmc-util.R")

model_data <- list(
  ti = 26,
  m = rbind(
    readRDS("rds/owls-example/capture-recapture-female-first-data.rds"),
    readRDS("rds/owls-example/capture-recapture-female-adult-data.rds")
  ),
  mM = rbind(
    readRDS("rds/owls-example/capture-recapture-male-first-data.rds"),
    readRDS("rds/owls-example/capture-recapture-male-adult-data.rds")
  )
)

model_data$r <- rowSums(model_data$m)
model_data$rM <- rowSums(model_data$mM)

n_chain <- 6
n_iter <- 2e4

capture_recapture_submodel <- jags.model(
  file = "scripts/owls-example/models/capture-recapture.bug",
  data = model_data,
  n.chains = n_chain,
  n.adapt = n_iter / 2
)

parameters <- c(
  "phij",
  "phia",
  "phijM",
  "phiaM",
  "p",
  "pM",  
  "v",
  "bp"
)

results <- coda.samples(
  model = capture_recapture_submodel,
  variable.names = parameters,
  n.iter = n_iter
)

# reflow into iters X chains X parameters
results_array <- mcmc_list_to_array(
  results
)

saveRDS(
  object = results_array,
  file = "rds/owls-example/capture-recapture-subposterior-samples.rds"
)
