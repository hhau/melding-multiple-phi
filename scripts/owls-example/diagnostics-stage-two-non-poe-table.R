source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

original_model_samples <- readRDS("rds/owls-example/original-ipm-samples.rds")
normal_approx_samples <- readRDS("rds/owls-example/melded-posterior-normal-approx-samples.rds")
melded_model_log_pooling_phi_samples <- readRDS("rds/owls-example/melded-phi-samples-log-pooling.rds")
melded_model_lin_pooling_phi_samples <- readRDS("rds/owls-example/melded-phi-samples-lin-pooling.rds")

# the original code makes a number of empty nodes through loops, so the static
# alpha_5 is not the alpha 5 we are interested in.
pars <- c(sprintf('v[%d]', 1 : 6), 'fec', sprintf('bp[%d]', 1 : 25))
vals <- c(sprintf('alpha[%d]', c(0 : 2, 4 : 6)), 'rho', sprintf('bp[%d]', 1 : 25))
names(vals) <- pars

orig_latex_names <- c(
  sprintf(r"{$\alpha_{%d}$}", c(0 : 2, 4 : 6)),
  r"{$\rho$}",
  sprintf(r"{$\alpha_{5, %d}$}", c(1 : 25))
)

write_diagnostics_to_file(
  samples_array = original_model_samples[, , pars],
  row_latex_names = orig_latex_names,
  output_filename = 'tex-input/owls-example/appendix-info/0011-orig-ipm-diagnostics.tex',
  caption = r'{Numerical diagnostics for the original IPM.}',
  label = 'owls-orig-ipm-diag'
)

normal_approx_pars <- c('fec', sprintf('v[%d]', c(1, 2, 6)))
normal_approx_latex_names <- c(
  r"{$\rho$}",
  sprintf(r"{$\alpha_{%d}$}", c(0, 2, 6))
)

write_diagnostics_to_file(
  samples_array = normal_approx_samples[, , normal_approx_pars],
  row_latex_names = normal_approx_latex_names,
  output_filename = 'tex-input/owls-example/appendix-info/0021-stage-two-normal-approx-diagnostics.tex',
  caption = r'{Numerical diagnostics for the normal approximation $\widehat{\pd}_{\text{meld}}$}',
  label = 'owls-stage-two-normal-approx-diag'
)

melded_latex_names <- c(
  sprintf(r"{$\alpha_{%d}$}", c(0, 2)),
  r"{$\rho$}"
)

write_diagnostics_to_file(
  samples_array = melded_model_log_pooling_phi_samples,
  row_latex_names = melded_latex_names,
  output_filename = 'tex-input/owls-example/appendix-info/0022-stage-two-log-pooling-diagnostics.tex',
  caption = r'{Numerical diagnostics for the melded posterior using logarithmic pooling $\pd_{\text{meld, log}}$}',
  label = 'owls-stage-two-log-diag'
)

write_diagnostics_to_file(
  samples_array = melded_model_lin_pooling_phi_samples,
  row_latex_names = melded_latex_names,
  output_filename = 'tex-input/owls-example/appendix-info/0023-stage-two-lin-pooling-diagnostics.tex',
  caption = r'{Numerical diagnostics for the melded posterior using linear pooling $\pd_{\text{meld, log}}$}',
  label = 'owls-stage-two-log-diag'
)
