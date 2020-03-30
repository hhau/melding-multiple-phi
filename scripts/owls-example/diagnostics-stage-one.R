source("scripts/common/plot-settings.R")

library(gt)
library(bayesplot)
library(rstan)
library(coda)
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
pars <- sprintf("v[%d]", 1 : 3)

# model 3
fecundity_diagnostics <- monitor(fecunditiy_submodel_samples)
fecundity_traceplot <- mcmc_trace(fecunditiy_submodel_samples)

ggsave_halfheight(
  filename = "plots/owls-example/stage-one-diagnostics-fecundity.pdf",
  plot = fecundity_traceplot
)

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

capture_recapture_traceplot <- mcmc_trace(
  capture_recapture_submodel_samples[, , rhat_max_index, drop = FALSE]
)

ggsave_halfheight(
  filename = "plots/owls-example/stage-one-diagnostics-capture-recapture.pdf",
  plot = capture_recapture_traceplot
)

phi_capture_recapture_diagnostics <- monitor(
  capture_recapture_submodel_samples[, , pars],
  print = FALSE
)

res <- dplyr::bind_rows(
  rownames_to_column(as.data.frame(phi_capture_recapture_diagnostics)),
  rownames_to_column(as.data.frame(fecundity_diagnostics))
) %>% 
  select(par = rowname, n_eff, Rhat, Bulk_ESS, Tail_ESS)

res <- kable(
  x = res,
  format = "latex",
  booktabs = TRUE
) %>%  
  kable_styling(latex_options = "striped") %>% 
  column_spec(1, "1cm")

cat(
  res,
  file = "tex-input/owls-example/appendix-info/0010-stage-one-diagnostics.tex" 
)