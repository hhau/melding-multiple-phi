library(rstan)

data <- readRDS(file = "rds/owls-example/fecundity-data.rds")

stan_data <- list(
  "T" = nrow(data),
  N_breeding_females = data$N_breeding_females,
  N_offspring = data$N_offspring
)

prefit <- stan_model(
  file = "scripts/owls-example/models/fecundity-model.stan"
)

model_fit <- sampling(
  object = prefit, 
  data = stan_data,
  chains = 6, 
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  warmup = 500,
  iter = 5e3 + 500,
  cores = 6
)

samples <- extract(model_fit, permuted = FALSE, pars = "rho")

saveRDS(
  object = samples,
  file = "rds/owls-example/fecundity-subposterior-samples.rds"
)
