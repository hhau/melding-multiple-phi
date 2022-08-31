library(rjags)
library(coda)
library(magrittr)

source("scripts/common/mcmc-util.R")

p1_samples <- readRDS("rds/owls-example/capture-recapture-subposterior-samples.rds")
p3_samples <- readRDS("rds/owls-example/fecundity-subposterior-samples.rds") %>% 
  as.numeric()

# figure out normal approximation covariances and means
phi1_samples <- p1_samples[, , sprintf("v[%d]", 1 : 2)] %>% 
  as.numeric() %>% 
  array(dim = c(prod(dim(p1_samples)[1 : 2]), 2))

p1_cov <- cov(phi1_samples)
p1_mean <- apply(phi1_samples, 2, mean)
p3_cov <- var(p3_samples)
p3_mean <- mean(p3_samples) 

n_approx_mean <- c(p1_mean, p3_mean)
n_approx_cov <- p1_cov %>% 
  cbind(rep(0, 2)) %>%
  rbind(c(rep(0, 2), p3_cov))

n_approx_prec <- solve(n_approx_cov)

model_data <- list(
  ti = 26,
  popcount = as.numeric(unlist(readRDS("rds/owls-example/count-data.rds"))),
  n_approx_mean = n_approx_mean,
  n_approx_prec = n_approx_prec
)

n_chain <- 6
n_iter <- 2e4

count_data_submodel <- jags.model(
  file = "scripts/owls-example/models/count-data-normal-approx.bug",
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
  file = "rds/owls-example/melded-posterior-normal-approx-samples.rds"
)
