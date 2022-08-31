library(GGally)
library(ggplot2)
library(parallel)
library(sn)
library(dplyr)
library(mvtnorm)
library(stringr)
library(future.apply)

plan(multisession, workers = 5)

source('scripts/common/plot-settings.R')
source('scripts/common/logger-setup.R')
source('scripts/common/setup-argparse.R')

parser$add_argument("--normal-approx-prior-estimates")
parser$add_argument("--raw-monte-carlo-prior-samples")
args <- parser$parse_args()

normal_prior_approx_parameters <- readRDS(args$normal_approx_prior_estimates)
prior_samples_raw <- readRDS(args$raw_monte_carlo_prior_samples)

n_icustays <- length(normal_prior_approx_parameters) # check the structure of either one of the sample files for this
n_prior_samples <- 20000

sample_from_normal_approx <- function(n_sims, pars) {
  event_indicator <- rbinom(n = n_sims, size = 1, prob = pars$mix_weight)
  n_ones <- sum(event_indicator == 1)
  n_cens <- sum(event_indicator == 0)

  ones_samples <- rmvnorm(
    n = n_ones,
    mean = pars$rf_event_mu,
    sigma = pars$rf_event_sigma_mat
  )

  colnames(ones_samples) <- c('event_time', 'k_i', 'eta_before', 'eta_after')

  if (n_cens > 0) {
    cens_samples <- rmvnorm(
      n = n_cens,
      mean = pars$censored_event_mu,
      sigma = pars$censored_event_sigma_mat
    )
  } else {
    cens_samples <- matrix(nrow = 0, ncol = length(pars$censored_event_mu))
  }

  colnames(cens_samples) <- c('k_i', 'eta_before', 'eta_after')

  ones_samples <- ones_samples %>%
      as_tibble() %>%
      mutate(event_indicator = 1)

  cens_samples <- cens_samples %>%
      as_tibble() %>%
      mutate(event_indicator = 0, event_time = NA)

  final_samples <- bind_rows(ones_samples, cens_samples) %>% mutate(
    event_time = brms::inv_logit_scaled(
      event_time,
      pars$lower_event_time_limit,
      pars$upper_event_time_limit
    ),
    k_i = brms::inv_logit_scaled(
      k_i,
      pars$lower_breakpoint_limit,
      pars$upper_breakpoint_limit
    ),
    eta_before = exp(eta_before),
    eta_after = exp(eta_after)
  ) %>%
    tidyr::replace_na(list(event_time = pars$upper_event_time_limit))
}

future_lapply(1 : n_icustays, future.seed = TRUE, function(icustay_index) {

  source('scripts/common/plot-settings.R')
  source('scripts/common/logger-setup.R')

  local_normal_approx_pars <- normal_prior_approx_parameters[[icustay_index]]
  samples_from_normal_approx <- sample_from_normal_approx(
    n_sims = n_prior_samples,
    pars = local_normal_approx_pars
  ) %>%
    mutate(type = 'normal approx')

  samples_from_monte_carlo <- prior_samples_raw[[icustay_index]] %>%
    select(event_time, k_i, eta_before, eta_after, event_indicator) %>%
    mutate(type = 'monte carlo')

  plot_tbl <- bind_rows(samples_from_normal_approx, samples_from_monte_carlo) %>%
    mutate(
      event_indicator = as.factor(event_indicator),
      type = as.factor(type),
      "italic(T)[i]" = event_time,
      "kappa[i]" = k_i,
      "eta[1 * ',' ~ i]^{b}" = eta_before,
      "eta[1 * ',' ~ i]^{a}" = eta_after,
      col_lab = paste0(event_indicator, '.', type)
    )

  p1 <- ggpairs(
    data = plot_tbl %>%
      filter(!is.na(event_indicator)),
    columns = 7 : 10,
    mapping = aes(colour = col_lab),
    lower = list(continuous = function(data, mapping ,...) {
      ggplot(mapping = mapping) +
        geom_point(
          data = data %>% filter(event_indicator == 0),
          pch = '.',
          alpha = 0.5
        ) +
        geom_density_2d(
          data = data %>% filter(event_indicator == 1),
          bins = 20,
          alpha = 0.5,
          n = 512 / 2
        ) +
        scale_colour_manual(
          values = c(
            '0.monte carlo' = blues[1],
            '0.normal approx' = blues[3],
            '1.monte carlo' = greens[1],
            '1.normal approx' = greens[4]
          ),
          labels = c(
            '0.monte carlo' = expression(italic('d')[italic(i)] == 0 * ',' ~ 'MC'),
            '0.normal approx' = expression(italic('d')[italic(i)] == 0 * ',' ~ 'N'),
            '1.monte carlo' = expression(italic('d')[italic(i)] == 1 * ',' ~ 'MC'),
            '1.normal approx' = expression(italic('d')[italic(i)] == 1 * ',' ~ 'N')
          )
        )
    }),
    diag = list(continuous = function(data, mapping, ...) {
      ggally_barDiag(data, mapping, bins = 25, alpha = 0.75 ,...) +
        scale_fill_manual(
          name = 'Sample type',
          values = c(
            '0.monte carlo' = blues[1],
            '0.normal approx' = blues[3],
            '1.monte carlo' = greens[1],
            '1.normal approx' = greens[4]
          ),
          labels = c(
            '0.monte carlo' = expression(italic('d')[italic(i)] == 0 * ',' ~ 'MC'),
            '0.normal approx' = expression(italic('d')[italic(i)] == 0 * ',' ~ 'NA'),
            '1.monte carlo' = expression(italic('d')[italic(i)] == 1 * ',' ~ 'MC'),
            '1.normal approx' = expression(italic('d')[italic(i)] == 1 * ',' ~ 'NA')
          )
        ) +
        theme(
         legend.text.align = 0
        )
    }),
    upper = list(continuous = function(data, mapping ,...) {
      ggplot(mapping = mapping) +
        geom_point(
          data = data %>% filter(event_indicator == 0),
          pch = '.',
          alpha = 0.5
        ) +
        geom_density_2d(
          data = data %>% filter(event_indicator == 1),
          bins = 20,
          alpha = 0.5,
          n = 512 / 2
        ) +
        scale_colour_manual(
          values = c(
            '0.monte carlo' = blues[1],
            '0.normal approx' = blues[3],
            '1.monte carlo' = greens[1],
            '1.normal approx' = greens[4]
          ),
          labels = c(
            '0.monte carlo' = expression(italic('d')[italic(i)] == 0 * ',' ~ 'MC'),
            '0.normal approx' = expression(italic('d')[italic(i)] == 0 * ',' ~ 'N'),
            '1.monte carlo' = expression(italic('d')[italic(i)] == 1 * ',' ~ 'MC'),
            '1.normal approx' = expression(italic('d')[italic(i)] == 1 * ',' ~ 'N')
          )
        ) +
        scale_x_log10() +
        scale_y_log10()
      }),
    title = bquote(italic('i')==.(icustay_index)),
    legend = c(1, 1),
    axisLabels = 'show',
    labeller = label_parsed
  )

  flog.info(
    msg = sprintf('Saving prior p2 plot: %d', icustay_index),
    name = base_filename
  )

  ggsave_fullpage(
    filename = sprintf("plots/mimic-example/p3-prior-pairs/pairs-%02d.png", icustay_index),
    plot = p1,
    adjust_height = -5,
  )
})
