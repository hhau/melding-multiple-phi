library(dplyr)
library(rstan)
library(kableExtra)

stage_two_samples <- readRDS("rds/owls-example/melded-posterior-samples.rds")
stage_two_diagnostics <-  monitor(stage_two_samples, print = FALSE)

# numerical diagnostics for the following parameters + the worst performing 
# other parameter
vars <- c("fec", "v[1]", "v[2]")

bulk_min_index <- which.min(stage_two_diagnostics$Bulk_ESS)
tail_min_index <- which.min(stage_two_diagnostics$Tail_ESS)
rhat_max_index <- which.max(stage_two_diagnostics$Rhat)

quantities_of_interest <- names(stage_two_samples[
  1,
  1,
  c(bulk_min_index, tail_min_index, rhat_max_index)
])

diag_vars <- unique(c(vars, quantities_of_interest))

# table output
stage_two_table <- as.data.frame(stage_two_diagnostics)[diag_vars, ] %>%
  tibble::rownames_to_column() %>%
  select(par = rowname, n_eff, Rhat, Bulk_ESS, Tail_ESS) %>%
  rename(
    "Parameter" = par,
    "$N_{\\text{eff}}$" = n_eff,
    "$\\widehat{R}$" = Rhat,
    "Bulk ESS" = Bulk_ESS,
    "Tail ESS" = Tail_ESS,
  ) 

parameter_recode_vector <- c(
  'v[1]' = '$\\alpha_{0}$',
  'v[2]' = '$\\alpha_{2}$',
  'v[6]' = '$\\alpha_{5}$',
  'rho' = '$\\rho$',
  'fec' = '$\\rho$'
) 

stage_two_table$Parameter <- stage_two_table$Parameter %>%
  recode(!!!parameter_recode_vector)

res <- kable(
    x = stage_two_table, 
    format = "latex",
    booktabs = TRUE,
    escape = FALSE
  ) %>%  
  kable_styling(latex_options = c("striped",  "hold_position")) %>% 
  column_spec(1, "2cm")

cat(
  res,
  file = "tex-input/owls-example/appendix-info/0020-stage-two-diagnostics.tex" 
)

