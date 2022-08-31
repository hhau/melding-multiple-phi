library(rstan)
library(kableExtra)
library(dplyr)
library(tibble)

# load model 1 and model 3 samples
capture_recapture_submodel_samples <- readRDS(
  "rds/owls-example/capture-recapture-subposterior-samples.rds"
)
fecunditiy_submodel_samples <- readRDS(
  "rds/owls-example/fecundity-subposterior-samples.rds"
)

# find then parameters within each model that have
# minimum ESS
# Maximum Rhat
# and plot the traces for each
pars <- sprintf("v[%d]", c(1, 2))

fecundity_diagnostics <- monitor(fecunditiy_submodel_samples, print = FALSE)

capture_recapture_diagnostics <- monitor(
  capture_recapture_submodel_samples[, , pars],
  print = FALSE
)

parameter_recode_vector <- c(
  'v[1]' = '$\\alpha_{0}$',
  'v[2]' = '$\\alpha_{2}$',
  'rho' = '$\\rho$'
)

res <- bind_rows(
  rownames_to_column(as.data.frame(capture_recapture_diagnostics)),
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
  digits = 2,
  booktabs = TRUE,
  escape = FALSE
) %>%
  kable_styling(latex_options = c("striped",  "hold_position")) %>%
  column_spec(1, "2cm")

cat(
  kable_res,
  file = "tex-input/owls-example/appendix-info/0010-stage-one-diagnostics.tex"
)
