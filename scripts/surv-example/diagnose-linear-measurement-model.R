library(bayesplot)
library(ggplot2)

# TODO: Decide what to include in the appendix
# keep the parcoord + pairs for now as they help understand what's actually
# wrong with poorly fitting models
 
source("scripts/common/plot-settings.R")

model_output <- readRDS(
  file = "rds/surv-example/linear-measurement-output.rds"
)

div_style <- parcoord_style_np(
  div_color = "green", 
  div_size = 0.5, 
  div_alpha = 0.9
)

mcmc_parcoord(
  model_output$samples,
  transformations = function(x) {(x - mean(x)) / sd(x)},
  size = 0.25,
  alpha = 0.05,
  np = model_output$nuts_params, 
  np_style = div_style
) +
  theme(
    axis.text.x = element_text(angle = 90) 
  )

mcmc_pairs(
  model_output$samples,
  pars = c("lp__", "sigma_beta_one", "sigma_beta_zero", "beta_one[8]"),
  np = model_output$nuts_params
)
