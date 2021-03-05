library(survival)
library(survminer)
library(dplyr)
library(tidyverse)
library(pbapply)
library(latex2exp)

source("scripts/common/plot-settings.R")
source("scripts/surv-example/GLOBALS.R")

submodel_2_data <- readRDS("rds/surv-example/submodel-two-simulated-data.rds")
psi_2_samples <- readRDS("rds/surv-example/stage-two-psi-2-samples.rds")
phi_23_samples <- readRDS("rds/surv-example/stage-two-phi-23-samples.rds")
phi_12_samples <- readRDS("rds/surv-example/stage-two-phi-12-samples.rds")
psi_2_fixed_samples <- readRDS("rds/surv-example/point-est-psi-2-samples.rds")
phi_23_point_est <- readRDS("rds/surv-example/stage-one-phi-23-posterior-median.rds")
psi_2_names <- c("theta_zero", "theta_one", "hazard_gamma", "alpha")

eta_zero_names <- grep("*,1]", names(phi_23_samples[1, 1,]), value = TRUE)
eta_one_names <- grep("*,2]", names(phi_23_samples[1, 1,]), value = TRUE)

n_plot <- 20
plot_t <- seq(from = 0.25, to = 1, length.out = n_plot)
n_iter <- dim(phi_12_samples)[1]
n_thin_per_chain <- 100
n_chain <- 5
thin_vec <- seq.int(from = 1 , to = n_iter, length.out = n_thin_per_chain) %>%
  round()

phi_12_thin_vec <- seq.int(1, n_iter, length.out = 5 * n_thin_per_chain) %>%
  round()

event_time_names <- names(phi_12_samples[1, 1,]) %>%
  grep("time", x = ., value = T)

event_indicator_names <- names(phi_12_samples[1, 1, ]) %>%
  grep("indicator", x = ., value = T)

# it would be nice to speed this up some, or at least cache the results
res <- pblapply(phi_12_thin_vec, cl = 6, function(iter_id) {
  df <- data.frame(
    event_time = phi_12_samples[iter_id, 1, event_time_names],
    event_indicator = phi_12_samples[iter_id, 1, event_indicator_names],
    iter_id = iter_id
  )
  surv_obj <- survfit(Surv(event_time, event_indicator) ~ 1, data = df) %>% 
    surv_summary() %>%
    as_tibble() %>%
    select(-c(std.err, upper, lower)) %>%
    mutate(iter_id = iter_id)
  
  base_df <- surv_obj[1, ] %>%
    mutate(time = 0, surv = 1, n.event = 0)
  
  min_t_df <- surv_obj[1, ] %>%
    mutate(time = 0.37, surv = 1, n.event = 0)
  
  res <- rbind(base_df, min_t_df, surv_obj)
}) 

plot_df <- bind_rows(res)
censor_df <- plot_df %>% 
  filter(n.censor > 0) %>%
  mutate(plot_alpha = (n.censor - min(n.censor)) / (diff(range(n.censor))))

base_plot <- ggplot(
  data = plot_df, 
  aes(x = time, y = surv, group = iter_id)
) +
  geom_step(alpha = 0.1) + 
  scale_x_continuous(limits = c(0.36, 1)) +
  geom_point(
    data = censor_df,
    inherit.aes = FALSE,
    mapping = aes(x = time, y = surv),
    alpha = censor_df %>% pull(plot_alpha),
    shape = "|",
    col = "black",
    size = 3
  )

log_surv_prob <- function(
  t,
  baseline_val,
  hazard_gamma,
  theta_zero,
  theta_one,
  alpha,
  eta_zero,
  eta_one
) {
  term_one <- -hazard_gamma * exp(
    theta_zero + baseline_val * theta_one + alpha * eta_zero
  )
  
  integrand_stable <- function(u) {
    sign <- (alpha / abs(alpha)) * (eta_one / abs(eta_one))
    lf <- log(abs(alpha)) + log(abs(eta_one)) - log(hazard_gamma) +  
      hazard_gamma * log(u) + alpha * eta_one * u
    res <- sign * exp(lf)
    return(res)
  }
  
  I <- integrate(
    integrand_stable, 
    lower = 0, 
    upper = t,
    rel.tol = 1e-10
  )$value
  
  W  <- exp(
    hazard_gamma * log(t) + 
    alpha * eta_one * t - 
    log(hazard_gamma)
  )
  res <- term_one * (W - I)
  return(res)
}

# test_list <- as.list(surv_prob_samples[1, 4 : 10]) %>%
#   c(list(t = plot_t[1]))
# res <- do.call(log_surv_prob, test_list)

psi_2_thinned_samples_tbl <- matrix(
  psi_2_samples[thin_vec, , ], 
  nrow = n_thin_per_chain * n_chain,
  ncol = dim(psi_2_samples)[3],
  dimnames = list(NULL, names(psi_2_samples[1, 1, ]))
) %>%
  as_tibble() %>%
  mutate(
    chain_id = rep(1 : n_chain, each = n_thin_per_chain),
    sample_id = rep(thin_vec, times = n_chain)
  )

phi_23_thinned_samples_tbl <- tibble(
  patient_id = rep(1 : n_patients, each = n_thin_per_chain * n_chain),
  chain_id = rep(rep(1 : n_chain, each = n_thin_per_chain), times = n_patients),
  sample_id = rep(thin_vec, times = n_chain * n_patients),
  eta_zero = phi_23_samples[thin_vec, , eta_zero_names] %>% as.numeric(),
  eta_one = phi_23_samples[thin_vec, , eta_one_names] %>% as.numeric()
)

surv_prob_samples <- phi_23_thinned_samples_tbl %>%
  left_join(psi_2_thinned_samples_tbl) %>%
  left_join(submodel_2_data)

haz_df <- pblapply(plot_t, cl = 6, function(a_t) {
  surv_prob_samples %>%
    rowwise() %>%
    mutate(
      plot_t = a_t,
      surv_prob = log_surv_prob(
        t = a_t,
        baseline_val = baseline_val,
        hazard_gamma = hazard_gamma,
        theta_zero = theta_zero,
        theta_one = theta_one,
        alpha = alpha,
        eta_zero = eta_zero,
        eta_one = eta_one
      ) %>% 
        exp()
    )
}) %>%
  bind_rows()

haz_plot_tbl <- haz_df %>%
  # filter(surv_prob > 0.01) %>% ## temporary hack because the numerics r wrong
  group_by(plot_t) %>%
  summarise(
    mean = mean(surv_prob),
    median = median(surv_prob),
    lower = quantile(surv_prob, 0.9),
    upper = quantile(surv_prob, 0.1)
  ) %>%
  mutate(method = "chained")

with_melded_interval <- base_plot +
  geom_line(
    data = haz_plot_tbl,
    inherit.aes = FALSE,
    aes(x = plot_t, y = mean), 
    col = blues[2]
  ) + 
  geom_ribbon(
    data = haz_plot_tbl,
    inherit.aes = FALSE,
    aes(x = plot_t, ymin = lower, ymax = upper), 
    alpha = 0.2, 
    col = blues[2]
  ) 

## now add the fully fixed version
psi_2_fixed_thinned_samples_tbl <- matrix(
  psi_2_fixed_samples[thin_vec, , psi_2_names], 
  nrow = n_thin_per_chain * n_chain,
  ncol = length(psi_2_names),
  dimnames = list(NULL, psi_2_names)
) %>%
  as_tibble() %>%
  mutate(
    chain_id = rep(1 : n_chain, each = n_thin_per_chain),
    sample_id = rep(thin_vec, times = n_chain)
  )

phi_23_fixed_thinned_samples_tbl <- tibble(
  patient_id = rep(1 : n_patients, each = n_thin_per_chain * n_chain),
  chain_id = rep(rep(1 : n_chain, each = n_thin_per_chain), times = n_patients),
  sample_id = rep(thin_vec, times = n_chain * n_patients),
  eta_zero = rep(phi_23_point_est$long_beta_zero, each = n_thin_per_chain * n_chain),
  eta_one = rep(phi_23_point_est$long_beta_one, each = n_thin_per_chain * n_chain)
)

surv_prob_fixed_samples <- phi_23_thinned_samples_tbl %>%
  left_join(psi_2_fixed_thinned_samples_tbl) %>%
  left_join(submodel_2_data)

haz_fixed_df <- pblapply(plot_t, cl = 6, function(a_t) {
  surv_prob_fixed_samples %>%
    rowwise() %>%
    mutate(
      plot_t = a_t,
      surv_prob = log_surv_prob(
        t = a_t,
        baseline_val = baseline_val,
        hazard_gamma = hazard_gamma,
        theta_zero = theta_zero,
        theta_one = theta_one,
        alpha = alpha,
        eta_zero = eta_zero,
        eta_one = eta_one
      ) %>% 
        exp()
    )
}) %>%
  bind_rows()

haz_fixed_plot_tbl <- haz_fixed_df %>%
  # filter(surv_prob > 0.01) %>% ## temporary hack because the numerics r wrong
  group_by(plot_t) %>%
  summarise(
    mean = mean(surv_prob),
    median = median(surv_prob),
    lower = quantile(surv_prob, 0.9),
    upper = quantile(surv_prob, 0.1)
  ) %>%
  mutate(method = "fixed")

both_haz_plot_tbl <- bind_rows(haz_plot_tbl, haz_fixed_plot_tbl)

p1 <- base_plot +
  geom_line(
    data = both_haz_plot_tbl,
    inherit.aes = FALSE,
    aes(x = plot_t, y = mean, col = method)
  ) + 
  geom_ribbon(
    data = both_haz_plot_tbl,
    inherit.aes = FALSE,
    aes(x = plot_t, ymin = lower, ymax = upper, col = method), 
    alpha = 0.2
  ) +
  scale_colour_manual(
    values = c(
      chained = highlight_col,
      fixed = blues[1]
    ),
    labels = list(
      chained = "Chained Melding",
      fixed = TeX("Fix $\\phi_{1 \\bigcap 2}$ and $\\phi_{2 \\bigcap 3}$")
    )
  ) +
  labs(colour = "Posterior type") +
  ylab(expression("S"(italic(t)))) +
  xlab(expression(italic(t))) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

ggsave_halfheight(
  filename = "plots/surv-example/kaplan-meier-pc.pdf",
  plot = p1
)
