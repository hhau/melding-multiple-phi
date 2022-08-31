library(bayesplot)
library(ggplot2)
 
source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

model_output <- readRDS(
  file = "rds/surv-example/submodel-one-output.rds"
)

model_samples <- model_output$samples

# need to cut the generated quantities
sub_vec <- !grepl(
  "(event_time|plot_mu|event_indicator|lp__)", 
  x = names(model_samples[1, 1, ])
)

p_1 <- plot_worst_pars(
  model_samples[, , sub_vec],
  facet_name_value_pairs = c(
    r"{beta_zero\[}" = r"{beta\[0 * ',' ~}",
    r"{beta_one\[}" = r"{beta\[1 * ',' ~}",
    "mu_beta_zero" = r"{mu[beta * ',' ~ 0]}",
    "mu_beta_zero" = r"{mu[beta * ',' ~ 1]}",
    "sigma_beta_zero" = r"{sigma[beta * ',' ~ 0]}",
    "sigma_beta_one" = r"{sigma[beta * ',' ~ 1]}",
    "sigma_y" = r"{sigma['y']}" 
  )
)

ggsave_halfheight(
  filename = "plots/surv-example/stage-one-submodel-one-diags.png",
  plot = p_1
)
