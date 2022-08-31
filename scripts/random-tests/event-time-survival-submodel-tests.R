library(tidyverse)
library(scales)

source('scripts/common/plot-settings.R')

generate_prior_sample <- function() {
  res <- list(
    hazard_gamma = abs(rnorm(n = 1)),
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

pointwise_unnormalised_target <- function(t, prior_sample, return_log = TRUE) {
  with(prior_sample, {
    log_hazard_term <- log(hazard_gamma) + ((hazard_gamma - 1) * log(t)) + (baseline_term)
    common_log_surv_term <- -exp(baseline_term)
    
    if (t < k_i) {
      longitudinal_term <- alpha * eta_before
      log_hazard_term <- log_hazard_term + longitudinal_term
      
      temp_term_1 <- exp(longitudinal_term) * t^(hazard_gamma)
      log_surv_term <- common_log_surv_term * temp_term_1
    } else {
      longitudinal_term_1 <- alpha * eta_before
      longitudianl_term_2 <- alpha * eta_after
      log_hazard_term <- log_hazard_term + longitudianl_term_2
      
      temp_term_2 <- exp(longitudinal_term_1) * (k_i^(hazard_gamma))
      temp_term_3 <- exp(longitudianl_term_2) * (t^(hazard_gamma) - k_i^(hazard_gamma))
      log_surv_term <- common_log_surv_term * (temp_term_2 + temp_term_3)
    }
    
    log_target <- log_hazard_term + log_surv_term
    
    if (return_log) {
      return(log_target)
    } else {
      return(exp(log_target))
    }
  })
}

get_normalising_constant <- function(
  prior_sample, 
  int_lower = 0, 
  int_upper = 100
) {
  temp_f <- function(x) {
    pointwise_unnormalised_target(x, prior_sample, return_log = FALSE)
  } 
    
  temp_f <- Vectorize(temp_f)
  normalising_constant <- integrate(
    f = temp_f,
    lower = int_lower,
    upper = int_upper,
    subdivisions = 1e3,
    stop.on.error = FALSE
  )$value
  
  return(normalising_constant)
}

plot_t <- exp(seq(from = log(0.00001), to = log(100), length.out = 1000))
plot_tbl <- lapply(1 : 24, function(iter) {
  a_sample <- generate_prior_sample()
  nc <- get_normalising_constant(a_sample)
  
  anon_f <- function(t) {
    pointwise_unnormalised_target(t, a_sample, return_log = FALSE)
  }
  
  res <- tibble(
    iter = iter,
    t = plot_t,
    q_val = sapply(plot_t, anon_f),
    p_val = q_val / nc,
    k_i = a_sample$k_i,
    hazard_gamma = a_sample$hazard_gamma,
    alpha = a_sample$alpha,
    baseline_term = a_sample$baseline_term
  )
  
  return(res)
}) %>% 
  bind_rows()

label_maker <- function(h, a, b) {
  sprintf(
    "gamma == %.3f * ',' ~  alpha == %.3g * ',' ~  bold(x)[i] * theta == %.3f",
    h,
    a,
    b
  )
}

p1 <- plot_tbl %>% 
  filter(!(is.na(q_val) | is.infinite(q_val))) %>% 
  mutate(
    plot_label = label_maker(hazard_gamma, alpha, baseline_term)
  ) %>% 
  ggplot() + 
  geom_line(aes(x = t, y = q_val), alpha = 0.5) +
  geom_vline(aes(xintercept = k_i), linetype = 'dotted', alpha = 0.8, col = highlight_col) +
  facet_wrap(vars(plot_label), scales = 'free_y', ncol = 3, labeller = label_parsed) +
  scale_y_continuous(labels = label_scientific(digits = 3)) +
  scale_x_log10() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(
      size = rel(0.8),
      margin = margin(2, 4, 2, 4, "pt")
    ),
    axis.text = element_text(size = rel(0.8))
  ) +
  xlab(expression('T'[i])) + 
  ylab(expression('Unnormalised target:' ~ 'q'('T'[i] ~ '|' ~ psi[2], ~ phi[23])))

ggsave_base(
  filename = 'plots/mimic-example/test-event-time-prior-target.pdf',
  plot = p1,
  height = display_settings$full_page_plot_height * (23 / 15),
  width = display_settings$full_page_plot_width * (23 / 15)
)
