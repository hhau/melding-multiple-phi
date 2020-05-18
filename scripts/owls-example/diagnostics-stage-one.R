# Thu 14 May 17:17:55 2020 - this file is poor because it does 2 things, and it
# shouldn't. Also the `monitor` call takes for ever, and should be saved out
# as a separate rds file.

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

library(bayesplot)
library(rstan)
library(kableExtra)
library(dplyr)
library(tibble)
library(parallel)

# load model 1 and model 3 samples 
capture_recapture_submodel_samples <- readRDS("rds/owls-example/capture-recapture-subposterior-samples.rds")
fecunditiy_submodel_samples <- readRDS("rds/owls-example/fecundity-subposterior-samples.rds")

# find then parameters within each model that have
# minimum ESS
# Maximum Rhat
# and plot the traces for each
pars <- sprintf("v[%d]", 1 : 2)

# model 3
fecundity_traceplot <- mcmc_combo(
  fecunditiy_submodel_samples,
  combo = c("trace", "rank_overlay")
)

ggsave_halfheight(
  filename = "plots/owls-example/stage-one-diagnostics-fecundity.png",
  plot = fecundity_traceplot
)

fecundity_diagnostics <- monitor(fecunditiy_submodel_samples)

capture_recapture_diagnostics <- mclapply(
  1 : dim(capture_recapture_submodel_samples)[3],
  mc.cores = 6,
  function(x) {
    monitor(capture_recapture_submodel_samples[, , x], print = FALSE) 
  }
) %>% 
  bind_rows() 

bulk_min_index <- which.min(capture_recapture_diagnostics$Bulk_ESS)
bulk_max_index <- which.max(capture_recapture_diagnostics$Bulk_ESS)
rhat_max_index <- which.max(capture_recapture_diagnostics$Rhat)

capture_recapture_parameter_names <- names(capture_recapture_submodel_samples[1, 1, ])
bulk_min_name <- capture_recapture_parameter_names[bulk_min_index]
rhat_max_name <- capture_recapture_parameter_names[rhat_max_index]

capture_recapture_traceplot <- mcmc_trace(
  capture_recapture_submodel_samples[, , c(pars, bulk_min_name, rhat_max_name), drop = FALSE]
)

ggsave_halfheight(
  filename = "plots/owls-example/stage-one-diagnostics-capture-recapture.png",
  plot = capture_recapture_traceplot
)

phi_capture_recapture_diagnostics <- monitor(
  capture_recapture_submodel_samples[, , pars],
  print = FALSE
)

parameter_recode_vector <- c(
  'v[1]' = '$\\alpha_{0}$',
  'v[2]' = '$\\alpha_{2}$',
  'rho' = '$\\rho$'
) 

res <- dplyr::bind_rows(
  rownames_to_column(as.data.frame(phi_capture_recapture_diagnostics)),
  rownames_to_column(as.data.frame(fecundity_diagnostics))
) %>% 
  select(par = rowname, n_eff, Rhat, Bulk_ESS, Tail_ESS) %>%
  rename(
    "Parameter" = par,
    "$N_{\\text{eff}}$" = n_eff,
    "$\\widehat{R}$" = Rhat,
    "Bulk ESS" = Bulk_ESS,
    "Tail ESS" = Tail_ESS,
  )

res$Parameter <- res$Parameter %>%
  recode(!!!parameter_recode_vector) 

kable_res <- kable(
  x = res,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE
) %>%  
  kable_styling(latex_options = c("striped",  "hold_position")) %>% 
  column_spec(1, "2cm")

cat(
  kable_res,
  file = "tex-input/owls-example/appendix-info/0010-stage-one-diagnostics.tex" 
)
