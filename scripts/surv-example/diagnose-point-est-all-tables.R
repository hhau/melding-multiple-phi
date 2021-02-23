source("scripts/common/mcmc-util.R")
source("scripts/surv-example/GLOBALS.R")

# psi_2_samples
point_est <- readRDS("rds/surv-example/point-est-psi-2-samples.rds")
point_est_1_meld_23 <- readRDS("rds/surv-example/point-est-1-meld-23-psi-2-samples.rds")
point_est_3_meld_12 <- readRDS("rds/surv-example/point-est-3-meld-12-psi-2-samples.rds")

named_vec <- c(
  "theta_zero" = r"{$\theta_{0}$}",
  'theta_one' = r"{$\theta_{1}$}",
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

# the phis
point_est_3_meld_12_phi_12_samples <- readRDS(
  "rds/surv-example/point-est-3-meld-12-phi-12-samples.rds"
)

point_est_1_meld_23_phi_23_samples <- readRDS(
  "rds/surv-example/point-est-1-meld-23-phi-23-samples.rds"
)

write_diagnostics_to_file(
  samples_array = point_est_3_meld_12_phi_12_samples,
  row_latex_names = c(
    sprintf(r"{$T_{%d}$}", 1 : n_patients),
    sprintf(r"{$\delta_{%d}$}", 1 : n_patients)
  ),
  output_filename = "tex-input/surv-example/0103-point-est-3-meld-12-phi-12-diag.tex",
  caption = r"{Numerical diagnostics for $\phi_{1 \cap 2}$ when $\phi_{2 \cap 3}$ is fixed and submodels $\pd_{1}$ and $\pd_{2}$ are melded.}"
)

write_diagnostics_to_file(
  samples_array = point_est_1_meld_23_phi_23_samples,
  row_latex_names = c(
    sprintf(r"{$\eta_{0, %d}$}", 1 : n_patients),
    sprintf(r"{$\eta_{1, %d}$}", 1 : n_patients)
  ),
  output_filename = "tex-input/surv-example/0104-point-est-1-meld-23-phi-23-diag.tex",
  caption = r"{Numerical diagnostics for $\phi_{2 \cap 3}$ when $\phi_{1 \cap 2}$ is fixed and submodels $\pd_{2}$ and $\pd_{3}$ are melded.}"
)
