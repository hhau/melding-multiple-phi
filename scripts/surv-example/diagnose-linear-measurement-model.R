library(bayesplot)
library(ggplot2)

# TODO: Decide what to include in the appendix
# keep the parcoord + pairs for now as they help understand what's actually
# wrong with poorly fitting models
 
source("scripts/common/plot-settings.R")

model_output <- readRDS(
  file = "rds/surv-example/submodel-one-output.rds"
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
  np_style = div_style,
  pars = vars(-contains("event"), -contains("plot"))
) +
  theme(
    axis.text.x = element_text(angle = 90) 
  )

pairs_plot_output <- mcmc_pairs(
  model_output$samples,
  pars = c(
    "lp__",
    "sigma_alpha_zero",
    "sigma_alpha_one",
    "sigma_beta_one",
    "sigma_noise_x",
    "sigma_noise_t",
    "sigma_total",
    "mu_alpha_zero",
    "mu_alpha_one",
    "mu_beta_one",
    "alpha_zero[7]",
    "alpha_one[7]",
    "beta_one[7]"
  ),
  transformations = list(
    sigma_alpha_zero = "log",
    sigma_alpha_one = "log",
    sigma_beta_one = "log",
    sigma_noise_x = "log",
    sigma_noise_t = "log",
    sigma_total = "log"
  ),
  np = model_output$nuts_params,
  off_diag_args = list(alpha = 0.3)
)

ggsave(
  filename = "plots/surv-example/reg-model-diagnostics.png",
  plot = pairs_plot_output,
  width = 45,
  heigh = 45,
  units = "cm"
)
