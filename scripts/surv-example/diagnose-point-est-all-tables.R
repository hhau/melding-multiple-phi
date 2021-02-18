source("scripts/common/mcmc-util.R")

point_est <- readRDS("rds/surv-example/point-est-psi-2-samples.rds")
point_est_1_meld_23 <- readRDS("rds/surv-example/point-est-1-meld-23-psi-2-samples.rds")
point_est_3_meld_12 <- readRDS("rds/surv-example/point-est-3-meld-12-psi-2-samples.rds")

named_vec <- c(
  "beta_zero" = r"{$\theta_{0}$}",
  'beta_one' = r"{$\theta_{1}$}",
  "hazard_gamma" = r"{$\gamma$}",
  "alpha" = r"{$\alpha$}"
)

par_names <- names(named_vec)

write_diagnostics_to_file(
  samples_array = point_est[, , par_names],
  row_latex_names = named_vec,
  output_filename = "tex-input/surv-example/0100-point-est-diag.tex",
  caption = r"{Numerical diagnostics for the posterior of $\psi_{2}$ obtained using point estimates for both $\phi_{1 \cap 2}$ and $\phi_{2 \cap 3}$ in the survival example.}",
)

write_diagnostics_to_file(
  samples_array = point_est_1_meld_23,
  row_latex_names = named_vec,
  output_filename = "tex-input/surv-example/0101-point-est-1-meld-23-diag.tex",
  caption = r"{Numerical diagnostics for the posterior of $\psi_{2}$ obtained using the point estimate for $\phi_{1 \cap 2}$ in the survival example.}",
)

write_diagnostics_to_file(
  samples_array = point_est_3_meld_12,
  row_latex_names = named_vec,
  output_filename = "tex-input/surv-example/0102-point-est-3-meld-12-diag.tex",
  caption = r"{Numerical diagnostics for the posterior of $\psi_{2}$ obtained using the point estimate for $\phi_{2 \cap 3}$ in the survival example.}",
)
