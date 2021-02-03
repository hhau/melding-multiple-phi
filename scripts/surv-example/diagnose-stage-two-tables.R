library(knitr)
library(kableExtra)
library(dplyr)

source("scripts/common/mcmc-util.R")

phi_12_samples <- readRDS("rds/surv-example/stage-two-phi-12-samples.rds")
phi_23_samples <- readRDS("rds/surv-example/stage-two-phi-23-samples.rds")
psi_2_samples <- readRDS("rds/surv-example/stage-two-psi-2-samples.rds")

write_diagnostics_to_file(
  samples_array = phi_12_samples,
  row_latex_names = sprintf(r"($t_{%d}$)", 1 : 16),
  output_filename = "tex-input/surv-example/0090-stage-two-phi-12-diag.tex",
  caption = r"(Stage two numeric diagnostics for $\phi_{1 \cap 2}$ in the survival example)",
  label = "surv-stage-two-diag-phi-12"
)

write_diagnostics_to_file(
  samples_array = phi_23_samples,
  row_latex_names = sprintf(r"($\beta_{0, %d}$)", 1 : 16),
  output_filename = "tex-input/surv-example/0091-stage-two-phi-23-diag.tex",
  caption = r"(Stage two numeric diagnostics for $\phi_{2 \cap 3}$ in the survival example)",
  label = "surv-stage-two-diag-phi-23"
)

write_diagnostics_to_file(
  samples_array = psi_2_samples,
  row_latex_names = c(
    r"($\beta_{0}$)", 
    r"($\beta_{1}$)", 
    r"($\gamma$)", 
    r"($\alpha$)"
  ),
  output_filename = "tex-input/surv-example/0092-stage-two-psi-2-diag.tex",
  caption = r"(Stage two numeric diagnostics for $\psi_{2}$ in the survival example)",
  label = "surv-stage-two-diag-psi-2"
)
