library(bayesplot)
library(ggplot2)
 
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

write_diagnostics_to_file(
  samples_array = model_samples[, , sub_vec],
  row_latex_names = c(
    sprintf(r"{$\beta_{0, %d}$}", 1 : 16),
    sprintf(r"{$\beta_{1, %d}$}", 1 : 16),
    r"{$\mu_{\beta, 0}$}",
    r"{$\mu_{\beta, 1}$}",
    r"{$\sigma_{\beta, 0}$}",
    r"{$\sigma_{\beta, 1}$}",
    r"{$\sigma_{y}$}"
  ),
  output_filename = "tex-input/surv-example/0080-submodel-one-numeric-diags.tex",
  caption = r"{Stage one numeric diagnostics for both $\psi_{1}$ and $\phi_{1 \cap 2}$ under submodel 1.}",
  label = "surv-stage-one-submodel-one"
)
