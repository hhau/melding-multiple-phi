# logging
# source("scripts/common/logger-setup.R")
source("scripts/common/mcmc-util.R")

# packages
library(coda)
library(abind)
library(nimble)
library(dplyr)

# flog.info("Loading data and samples", name = base_filename)
n_chains <- 6

# Read in model two data   
count_model_data <- readRDS("rds/owls-example/count-data.rds")

# Set up top level parameters
n_iter <- 2e4

# nimble setup
count_model_consts <- list(
  ti = 26
)

count_model_data <- list(
  popcount = as.numeric(unlist(readRDS("rds/owls-example/count-data.rds")))
)
# setup model as nimblecode
# consult models/count-data.bug for commented version of this model.
# have had to use the nimble form of notation for the truncation of the
# prior on v[ii]
count_model_code <- nimbleCode({
  logit(phij) <- v[1]          
  logit(phia) <- v[1] + v[2]
  log(im) <- v[6]

  for (ii in 1:6) {
    v[ii] ~ T(dnorm(0, sd = 100), -10, 10)
  } 

  fec ~ dunif(0, 10)
  
  for (ii in 1:50) {
    flat_p[ii] <- 1 / 50
  } 
  
  N1[1] ~ dcat(flat_p[1:50])
  NadSurv[1] ~ dcat(flat_p[1:50])
  Nadimm[1] ~ dcat(flat_p[1:50])

  for (tt in 2:ti) {
    mean1[tt - 1] <- 0.5 * fec * phij * Ntot[tt - 1]
    N1[tt] ~ dpois(mean1[tt - 1])
    mpo[tt - 1] <- Ntot[tt - 1] * im
    NadSurv[tt] ~ dbin(phia, Ntot[tt - 1])
    Nadimm[tt] ~ dpois(mpo[tt - 1])
  }

  for (tt in 1:ti) {
    Ntot[tt] <- NadSurv[tt] + Nadimm[tt] + N1[tt]
    popcount[tt] ~ dpois(Ntot[tt])
  }
})

source("scripts/owls-example/mcmc-nimble-functions.R")

count_model <- nimbleModel( 
  code = count_model_code, 
  name = "count_model",
  constants = count_model_consts,
  data = count_model_data,
  inits = count_model_inits_function()
)

count_model_initial_compile <- compileNimble(
  count_model
)

# make original mcmc_config 
count_model_config <- configureMCMC(
  count_model,
  print = TRUE
)
count_model_config$removeSampler(
  "v[1:2]"
)
count_model_config$removeSampler(
  "fec"
)
count_model_config$removeSampler(
  "N1[2:26]"
)
count_model_config$removeSampler(
  "NadSurv[2:26]"
)
count_model_config$removeSampler(
  "Nadimm[2:26]"
)

all_discrete_slice_nodes <- c(
  count_model$expandNodeNames('N1[2:26]'),
  count_model$expandNodeNames('NadSurv[2:26]'),
  count_model$expandNodeNames('Nadimm[2:26]')
)

for (node in all_discrete_slice_nodes) {
  count_model_config$addSampler(
    target = node,
    type = "positive_discrete_slice_sampler",
    control = list(
      adaptInterval = 400,
      maxContractions = 1000
    )
  )  
}

count_model_config$addSampler(
  target = "v[1:2]",
  type = "phi_one_two_proposal",
  control = list(capture_recapture_proposal = capture_recapture_proposal)
)
count_model_config$addSampler(
  target = "fec",
  type = "phi_two_three_proposal",
  control = list(fecunditiy_proposal = fecunditiy_proposal)
)

count_model_mcmc <- buildMCMC(
  count_model_config,
  project = count_model
)

compiled_count_model <- compileNimble(
  count_model_mcmc,
  project = count_model,
  showCompilerOutput = TRUE
)

mcmc_output <- runMCMC(
  compiled_count_model,
  niter = n_iter + 500,
  nburnin = 500,
  nchains = 6,
  inits = count_model_inits_function,
  samplesAsCodaMCMC = TRUE
)

# logfile_line_delineator()

saveRDS(
  file = "rds/owls-example/melded-posterior-samples.rds",
  object = mcmc_list_to_array(mcmc_output)
)

