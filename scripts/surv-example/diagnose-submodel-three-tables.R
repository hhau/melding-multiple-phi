source("scripts/surv-example/GLOBALS.R")
source("scripts/common/mcmc-util.R")

model_output <- readRDS(file = "rds/surv-example/submodel-three-output.rds")
model_samples <- model_output$samples

# need to cut the generated quantities
sub_vec <- !grepl("(plot_mu|lp__)", x = names(model_samples[1, 1, ]))

write_diagnostics_to_file(
  samples_array = model_samples[, , sub_vec],
  row_latex_names = c(
    sprintf(r"{$\beta_{0, %d}$}", 1 : n_patients),
    # sprintf(r"{$\beta_{1, %d}$}", 1 : 16),
    r"{$\mu_{\beta, 0}$}",
    # r"{$\mu_{\beta, 1}$}",
    r"{$\sigma_{\beta, 0}$}",
    # r"{$\sigma_{\beta, 1}$}",
    r"{$\sigma_{y}$}"
  ),
  output_filename = "tex-input/surv-example/0081-submodel-three-numeric-diags.tex",
  caption = r"{Stage one numeric diagnostics for both $\psi_{3}$ and $\phi_{2 \cap 3}$ under submodel 3, in the survival example.}",
  label = "surv-stage-one-submodel-three"
)
