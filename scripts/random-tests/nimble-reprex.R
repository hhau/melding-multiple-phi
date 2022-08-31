library(nimble)

# samples - to be used as a proposal
proposal_samples <- array(
  data = rnorm(n = 100 * 2),
  dim = c(100, 2),
  dimnames = list(
    NULL,
    c("mu", "theta")
  )
)

# constants and data
model_consts <- list(
  n_data = 15
)
model_data <- list(
  y = rnorm(n = model_consts$n_data, mean = 1, sd = 2)
)

# model code
model_code <- nimbleCode({
  mu ~ dnorm(0, sd = 2)
  sigma_y ~ T(dnorm(0, sd = 1), 0, )

  for (ii in 1 : n_data) {
    y[ii] ~ dnorm(mu, sd = sigma_y)
  }
})

model_inits <- list(
  mu = rnorm(n = 1),
  sigma_y = abs(rnorm(n = 1))
)

# proposal code
mu_proposer <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    internal_proposal_samples <- control$proposal_samples # maybe passing it in here helps?
  },
  run = function() {
    model_log_prob_initial <- getLogProb(model, calcNodes)
    
    proposal_index <- rcat(
      n = 1,
      prob = nimRep(1, each = 100)
    )
    proposal <- internal_proposal_samples[
      proposal_index, 
      'mu' # I suspect this, change to numeric if doesn't work
    ]
    model[[target]] <<- proposal
    model_log_prob_proposed <- calculate(model, calcNodes)
    
    log_MH_ratio <- model_log_prob_proposed - model_log_prob_initial

    u <- runif(1, 0, 1)
    if (u < exp(log_MH_ratio)) {
      jump <- TRUE
    } else { 
      jump <- FALSE
    }
    if (jump) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list(
    reset = function() {}
  )
)

# compilation step
model <- nimbleModel(
  code = model_code,
  name = "model",
  constants = model_consts,
  data = model_data,
  inits = model_inits
)
compiled_model <- compileNimble(
  model
)

# mcmc_config setup 
model_config <- configureMCMC(
  model
)
model_config$removeSampler(
  "mu"
)
model_config$addSampler(
  target = "mu",
  type = "mu_proposer",
  control = list(
    proposal_samples = proposal_samples
  )
)
new_mcmc <- buildMCMC(
  model_config,
  project = model
)

# error generating step
compiled_model_new_mcmc <- compileNimble(
  new_mcmc,
  project = model,
  showCompilerOutput = TRUE
)

# run the thing
compiled_model_new_mcmc$run(1e3)
samples <- as.matrix(compiled_model_new_mcmc$mvSamples)

# reweight the initial
output_index_vector <- match(
  x = samples[, "mu"],
  table = proposal_samples[, "mu"]
)

resulting_posterior <- proposal_samples[
  output_index_vector,
]

