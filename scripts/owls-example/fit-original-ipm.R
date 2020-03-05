# library(nimble)
library(rjags)
library(coda)

# read in the data
model_data <- list(
  time = 1 : 25,
  ti = 26,
  nestlings = readRDS("rds/owls-example/fecundity-data.rds")[1 : 25,'N_offspring'],
  sample.size = readRDS("rds/owls-example/fecundity-data.rds")[1 : 25, 'N_breeding_females'],
  popcount = as.numeric(unlist(readRDS("rds/owls-example/count-data.rds"))),
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
# fit th1e model
original_model <- jags.model(
  file = "scripts/owls-example/models/original-ipm.bug",
  data = model_data,
  n.chains = n_chain,
  n.adapt = n_iter / 2
)

parameters <- c(
  "phij",
  "phia",
  "phijM",
  "phiaM",
  "fec",
  "im",
  "lambda",
  "p",
  "pM",
  "NadSurv",
  "Nadimm",
  "Ntot",
  "MEPHJUF",
  "MEPHADF",
  "MEFE",
  "MEIM_H",
  "MEIM_L",
  "MEPHJUM",
  "MEPHADM",
  "v",
  "bp"
)

results <- coda.samples(
  model = original_model,
  variable.names = parameters,
  n.iter = n_iter
)

# reflow into iters X chains X parameters
n_par <- ncol(results[[1]])
par_names <- dimnames(results[[1]])[[2]]

results_array <- array(
  data = NA,
  dim = c(n_iter, n_chain, n_par),
  dimnames = list(
    NULL,
    sprintf("chain_%d", 1 : n_chain),
    par_names
  )
)

for (chain_index in 1 : n_chain) {
  results_array[, chain_index, ] <- results[[chain_index]]
}

saveRDS(
  object = results_array,
  file = "rds/owls-example/original-ipm-results.rds"
)
