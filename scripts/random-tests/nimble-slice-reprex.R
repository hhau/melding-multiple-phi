library(nimble)

# model code
model_code <- nimbleCode({
  x ~ dpois(1.0)
})

# inits
model_inits <- list(
  x = 1
)

# initial setup and compile
model <- nimbleModel(
  code = model_code,
  name = "model",
  inits = model_inits
)
compiled_model <- compileNimble(
  model
)

# mcmc config - delete - add
model_config <- configureMCMC(
  model
)
model_config$removeSampler(
  "x"
)
model_config$addSampler(
  target = "x",
  type = "slice"
)

# recompile
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

# run
compiled_model_new_mcmc$run(1e5)
