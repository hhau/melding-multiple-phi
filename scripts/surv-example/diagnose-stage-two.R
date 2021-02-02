library(bayesplot)
library(ggplot2)
library(knitr)
library(kableExtra)
library(dplyr)

source("scripts/common/plot-settings.R")

# read in psi_2, phi_12 and phi_23 samples
# for each parameter in those parameters I would
# like to plot the trace, rank, acf, and local ESS 

phi_12_samples <- readRDS("rds/surv-example/stage-two-phi-12-samples.rds")
phi_23_samples <- readRDS("rds/surv-example/stage-two-phi-23-samples.rds")
psi_2_samples <- readRDS("rds/surv-example/stage-two-psi-2-samples.rds")

p_1 <- mcmc_trace(phi_12_samples)
p_2 <- mcmc_trace(phi_23_samples)
p_3 <- mcmc_trace(psi_2_samples)

# need two functions that will
# 1 - make a 2x2 plot of diagnostic plots and save them with an appropriate
# name

# 2 - make an appropriate latex table of numerical diagnostics.
write_diagnostics_to_file <- function(
  samples_array,
  row_latex_names,
  output_filename,
  monitor_cols = c("mean", "sd", "Rhat", "Bulk_ESS", "Tail_ESS"),
  col_latex_names = c(
    "Mean",
    "Std dev", 
    r"($\widehat{R}$)", 
    r"($\widehat{\text{ESS}}_{B}$)", 
    r"($\widehat{\text{ESS}}_{T}$)"
  ),
  n_warmup = 1,
  n_digits = 3,
  ...
) {
  table <- rstan::monitor(samples_array, n_warmup, print = FALSE) %>%
    # as.data.frame() %>%
    select(!!!monitor_cols)
  
  rownames(table) <- row_latex_names
  formatted_kable <- kable(
    table, 
    format = "latex", 
    digits = 3, 
    col.names = col_latex_names,
    escape = FALSE,
    booktabs = T,
    linesep = "",
    ...
  ) %>%
    kable_styling(latex_options = "striped")
  
  invisible(cat(formatted_kable, file = output_filename))
}

write_diagnostics_to_file(
  samples_array = phi_12_samples,
  row_latex_names = sprintf(r"($t_{%d}$)", 1 : 16),
  output_filename = "tex-input/surv-example/0090-stage-two-phi-12-diag.tex",
  caption = r"(Stage two numeric diagnostics for $\phi_{1 \cap 2}$)",
  label = "surv-stage-two-diag-phi-12"
)

write_diagnostics_to_file(
  samples_array = phi_23_samples,
  row_latex_names = sprintf(r"($\beta_{0, %d}$)", 1 : 16),
  output_filename = "tex-input/surv-example/0091-stage-two-phi-23-diag.tex",
  caption = r"(Stage two numeric diagnostics for $\phi_{2 \cap 3}$)",
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
  caption = r"(Stage two numeric diagnostics for $\psi_{2}$)",
  label = "surv-stage-two-diag-psi-2"
)
