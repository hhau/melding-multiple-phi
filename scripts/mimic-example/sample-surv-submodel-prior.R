library(simsurv)
library(tibble)
library(dplyr)
library(magrittr)
library(GGally)
library(sn)
library(parallel)
library(rstan)
library(stringr)
library(mvtnorm)

source('scripts/common/plot-settings.R')
source('scripts/common/logger-setup.R')
source('scripts/common/setup-argparse.R')

parser$add_argument("--fluid-stan-data")
parser$add_argument("--pf-list-data")
parser$add_argument("--baseline-covariate-data")
parser$add_argument("--submodel-one-median-both")
parser$add_argument("--surv-prior-optim-stan-model")
args <- parser$parse_args()

fluid_stan_data <- readRDS(args$fluid_stan_data)
pf_list_data <- readRDS(args$pf_list_data)
baseline_covariate_data <- readRDS(args$baseline_covariate_data)
submodel_one_median_both <- readRDS(args$submodel_one_median_both)

prefit <- stan_model(args$surv_prior_optim_stan_model)

n_prior_samples <- 20000
n_icustays <- length(pf_list_data)
log_crude_event_rate <- submodel_one_median_both[[2]] %>% pull()

center_keep_intercept <- function(x) {
  p <- ncol(x)
  res <- scale(x[, 2 : p], center = TRUE, scale = TRUE)
  cbind(1, res)
}

baseline_data_x_mat <- model.matrix(
  ~ 1 + aniongap_median + bicarbonate_median + creatinine_median + chloride_median +
    glucose_median + hematocrit_median + hemoglobin_median + platelet_median +
    potassium_median + ptt_median + inr_median + pt_median + sodium_median +
    bun_median + wbc_median + gender + age_at_icu_adm,
  data = baseline_covariate_data %>% as.data.frame()
) %>%
  center_keep_intercept()

pointwise_unnormalised_log_hazard <- function(t, x, betas) {
  with(betas, {
    log_hazard_term <- log(hazard_gamma) + ((hazard_gamma - 1) * log(t)) + (baseline_term)

    if (t < k_i) {
      longitudinal_term <- alpha * eta_before
      log_hazard_term <- log_hazard_term + longitudinal_term
    } else {
      longitudianl_term_2 <- alpha * eta_after
      log_hazard_term <- log_hazard_term + longitudianl_term_2
    }

    log_dd_term <- log(dd_gamma)
    log_res <- matrixStats::logSumExp(c(log_hazard_term, log_dd_term))

    return(log_res)
  })
}

prior_samples_raw <- mclapply(1 : n_icustays, mc.cores = 5, function(icustay_index) {
  flog.info(
    sprintf(
      "mimic-ex: simulating prior samples from survival model for id: %d",
      icustay_index
    ),
    name = base_filename
  )

  indiv_list_data  <- pf_list_data[[icustay_index]]
  lower_event_time_limit <- indiv_list_data$boundary_knots[1]
  upper_event_time_limit <- indiv_list_data$boundary_knots[2]
  lower_breakpoint_limit <- fluid_stan_data$breakpoint_lower[icustay_index]
  upper_breakpoint_limit <- fluid_stan_data$breakpoint_upper[icustay_index]
  breakpoint_scale <- (upper_breakpoint_limit - lower_breakpoint_limit)

  generate_prior_sample <- function() {
    res <- list(
      hazard_gamma = rgamma(n = 1, shape = 9.05, rate = 8.72),
      w_i = baseline_data_x_mat[icustay_index, ],
      theta = c(
        rnorm(n = 1, mean = log_crude_event_rate, sd = 0.5),
        sn::rsn(n = 17, xi = 0, omega = 0.5, alpha = -1)
      ),
      alpha = sn::rsn(n = 1, xi = 0, omega = 0.5, alpha = -2),
      k_i = (rbeta(n = 1, 5.0, 5.0) * breakpoint_scale) + lower_breakpoint_limit,
      eta_before = rgamma(n = 1, shape = 1.53, rate = 0.24),
      eta_after = rgamma(n = 1, shape = 1.53, rate = 0.24),
      dd_gamma = rgamma(n = 1, shape = 2, rate = 2)
    )

    res$baseline_term <- (res$w_i %*% res$theta) %>% as.numeric()
    return(res)
  }

  many_df <- lapply(1 : n_prior_samples, function(x) {
    generate_prior_sample() %>%
      inset2('w_i', NULL) %>%
      inset2('theta', NULL) %>%
      as_tibble()
  })

  event_time_res <- lapply(many_df, function(z) {
    tryCatch(
      expr = {
        simsurv(
          loghazard = pointwise_unnormalised_log_hazard,
          x = tibble(id = icustay_index),
          betas = z,
          interval = c(1e-8, 500),
          maxt = upper_event_time_limit,
          rootsolver = 'uniroot',
          seed = 12931
        )
      },
      error = function(e) {
        e
      }
    )
  })

  samples_with_errors <- lapply(1 : n_prior_samples, function(sample_id) {
    if ('error' %in% class(event_time_res[[sample_id]])) {
      event_time <- NA
      event_indicator <- NA
      simsurv_status <- event_time_res[[sample_id]][['message']]
    } else {
      event_time <- event_time_res[[sample_id]][['eventtime']]
      simsurv_status <- 'Ok'

      if (event_time_res[[sample_id]][['status']] == 0) {
       event_indicator <- 2
      } else {
        event_indicator <- 1
      }
    }

    inner_res <- many_df[[sample_id]]
    inner_res$event_time <- event_time
    inner_res$event_indicator <- event_indicator
    inner_res$simsurv_status <- simsurv_status

    return(inner_res)
  }) %>%
    bind_rows() %>%
    mutate(
      sample_id = 1 : n(),
      icustay_index = icustay_index
    )
})

normal_prior_approx_parameters <- mclapply(1 : n_icustays, mc.cores = 5, function(icustay_index) {
  indiv_list_data  <- pf_list_data[[icustay_index]]
  surv_prior_samples <- prior_samples_raw[[icustay_index]] %>%
    filter(!is.na(event_time))

  stan_data <- list(
    n_prior_samples = nrow(surv_prior_samples),
    lower_event_time_limit = 0,
    upper_event_time_limit = indiv_list_data$boundary_knots[2],
    lower_breakpoint_limit = fluid_stan_data$breakpoint_lower[icustay_index],
    upper_breakpoint_limit = fluid_stan_data$breakpoint_upper[icustay_index],
    event_time_samples = surv_prior_samples$event_time,
    event_indicator_samples = surv_prior_samples$event_indicator,
    breakpoint_samples = surv_prior_samples$k_i,
    eta_before_samples = surv_prior_samples$eta_before,
    eta_after_samples = surv_prior_samples$eta_after
  )

  keep_trying <- TRUE
  attempts <- 1

  while (keep_trying) {
    optm_res <- optimizing(
      prefit,
      data = stan_data
    )

    if (optm_res$return_code == 0) {
      keep_trying <- FALSE
    }

    if (attempts %% 10 == 0) {
      flog.info(
        "icustay_index %d is on attempt %d -- inspect fit if possible",
        icustay_index,
        attempts
      )
    }

    attempts <- attempts + 1
  }

  par_values <- list(
    icustay_index = icustay_index,
    optimizing_return_code = optm_res$return_code,
    lower_event_time_limit = 0,
    upper_event_time_limit = indiv_list_data$boundary_knots[2],
    lower_breakpoint_limit = fluid_stan_data$breakpoint_lower[icustay_index],
    upper_breakpoint_limit = fluid_stan_data$breakpoint_upper[icustay_index],
    mix_weight = optm_res$par['mix_weight'],
    rf_event_mu = optm_res$par %>%
      names() %>%
      str_detect("rf_event_mu") %>%
      magrittr::extract(optm_res$par, .) %>%
      as.numeric(),
    rf_event_sigma_mat = optm_res$par %>%
      names() %>%
      str_detect("rf_event_sigma_mat") %>%
      magrittr::extract(optm_res$par, .) %>%
      matrix(data = ., nrow = 4, ncol = 4),
    dd_event_mu = optm_res$par %>%
      names() %>%
      str_detect("dd_event_mu") %>%
      magrittr::extract(optm_res$par, .) %>%
      as.numeric(),
    dd_event_sigma_mat = optm_res$par %>%
      names() %>%
      str_detect("dd_event_sigma_mat") %>%
      magrittr::extract(optm_res$par, .) %>%
      matrix(data = ., nrow = 3, ncol = 3)
  )
})

sample_from_normal_approx <- function(n_sims, pars) {
  event_indicator <- rbinom(n = n_sims, size = 1, prob = 1 - pars$mix_weight) + 1
  n_ones <- sum(event_indicator == 1)
  n_twos <- sum(event_indicator == 2)

  ones_samples <- rmvnorm(
    n = n_ones,
    mean = pars$rf_event_mu,
    sigma = pars$rf_event_sigma_mat
  )

  colnames(ones_samples) <- c('event_time', 'k_i', 'eta_before', 'eta_after')

  if (n_twos > 0) {
    twos_samples <- rmvnorm(
      n = n_twos,
      mean = pars$dd_event_mu,
      sigma = pars$dd_event_sigma_mat
    )
  } else {
    twos_samples <- matrix(nrow = 0, ncol = length(pars$dd_event_mu))
  }

  colnames(twos_samples) <- c('k_i', 'eta_before', 'eta_after')

  ones_samples <- ones_samples %>%
      as_tibble() %>%
      mutate(event_indicator = 1)

  twos_samples <- twos_samples %>%
      as_tibble() %>%
      mutate(event_indicator = 2, event_time = NA)

  final_samples <- bind_rows(ones_samples, twos_samples) %>% mutate(
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

# plot the relevant cols and samples by method
lapply(1 : n_icustays, function(icustay_index) {
  local_normal_approx_pars <- normal_prior_approx_parameters[[icustay_index]]
  samples_from_normal_approx <- sample_from_normal_approx(
    n_prior_samples,
    local_normal_approx_pars
  ) %>%
    mutate(type = 'normal approx')

  samples_from_monte_carlo <- prior_samples_raw[[icustay_index]] %>%
    select(event_time, k_i, eta_before, eta_after, event_indicator) %>%
    mutate(type = 'monte carlo')

  plot_tbl <- bind_rows(samples_from_normal_approx, samples_from_monte_carlo) %>%
    mutate(
      event_indicator = as.factor(event_indicator),
      type = as.factor(type)
    )

  p1 <- ggpairs(
    data = plot_tbl,
    columns = 1 : 4,
    mapping = aes(colour = interaction(event_indicator, type)),
    lower = list(continuous = function(data, mapping ,...) {
      ggally_points(data, mapping, alpha = 0.5, shape = ".", ...) +
        scale_x_log10() +
        scale_y_log10()
    }),
    diag = list(continuous = function(data, mapping, ...) {
      ggally_barDiag(data, mapping, bins = 25, alpha = 0.75 ,...) +
        scale_x_log10()
    }),
    upper = list(continuous = 'blank'),
    title = bquote(italic('i')==.(icustay_index)),
    legend = c(1, 1)
  )

  ggsave_base(
    filename = sprintf("plots/mimic-example/temp-pairs/pairs-%02d.png", icustay_index),
    plot = p1,
    height = 25,
    width = 25
  )
})

saveRDS(
  object = normal_prior_approx_parameters,
  file = args$output
)
