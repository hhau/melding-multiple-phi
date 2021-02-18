library(bayesplot)
library(ggplot2)
 
source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")
source("scripts/surv-example/GLOBALS.R")

model_output <- readRDS(
  file = "rds/surv-example/submodel-three-output.rds"
)

model_samples <- model_output$samples

# need to cut the generated quantities
sub_vec <- !grepl(
  "(plot_mu|lp__)", 
  x = names(model_samples[1, 1, ])
)

names_vec <- c(
  sprintf("mu_beta[%d]", 1 : n_long_beta),
  sprintf("sigma_beta[%d]", 1 : n_long_beta),
  sprintf("beta[%d,1]", 1 : n_patients),
  sprintf("beta[%d,2]", 1 : n_patients),
  "sigma_y"
)

values_vec <- c(
  sprintf(r"{mu[eta * ',' ~ %d]}", (1 : n_long_beta) - 1),
  sprintf(r"{sigma[eta * ',' ~ %d]}", (1 : n_long_beta) - 1),
  sprintf(r"{eta[%d * ',' ~ 0]}", 1 : n_patients),
  sprintf(r"{eta[%d * ',' ~ 1]}", 1 : n_patients),
  r"{sigma['y']}"
)

names(values_vec) <- names_vec

rep_func <- function(x) {
  values_vec[x]
}

p_1 <- plot_worst_pars(
  model_samples[, , sub_vec],
  facet_name_value_pairs = rep_func
)

ggsave_halfheight(
  filename = "plots/surv-example/stage-one-submodel-three-diags.png",
  plot = p_1
)
