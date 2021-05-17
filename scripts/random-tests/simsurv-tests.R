library(simsurv)
library(tibble)
library(dplyr)

generate_prior_sample <- function() {
  res <- list(
    hazard_gamma = rgamma(n = 1, shape = 2, rate = 1),
    w_i = rnorm(n = 19),
    theta = rnorm(n = 19),
    alpha = rnorm(n = 1) / 2000 ,
    k_i = abs(rnorm(n = 1, mean = 4.5)),
    eta_before = abs(rnorm(n = 1, mean = 5e3, sd = 1e3)),
    eta_after = abs(rnorm(n = 1, mean = 5e3, sd = 1e3))
  )
  
  res$baseline_term <- as.numeric(res$w_i %*% res$theta)
  return(res)
}

pointwise_unnormalised_log_hazard <- function(t, x, betas) {
  with(betas, {
    log_hazard_term <- log(hazard_gamma) + ((hazard_gamma - 1) * log(t)) + (baseline_term)
    common_log_surv_term <- -exp(baseline_term)
    
    if (t < k_i) {
      longitudinal_term <- alpha * eta_before
      log_hazard_term <- log_hazard_term + longitudinal_term
    } else {
      longitudianl_term_2 <- alpha * eta_after
      log_hazard_term <- log_hazard_term + longitudianl_term_2
    }
    return(log_hazard_term)
  })
}

n_sims <- 200

many_df <- lapply(1 : n_sims, function(x) {
  generate_prior_sample() %>% 
    inset2('w_i', NULL) %>% 
    inset2('theta', NULL) %>% 
    as_tibble()
}) %>% 
  bind_rows()

simsurv(
  loghazard = pointwise_unnormalised_log_hazard,
  x = data.frame(ids = 1 : n_sims),
  betas = many_df,
  interval = c(1e-8, 2e7)
)
