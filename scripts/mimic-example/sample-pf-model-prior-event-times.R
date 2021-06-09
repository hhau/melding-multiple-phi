library(dplyr)
library(splines2)
library(rootSolve)
library(abind)
library(parallel)
library(magrittr)
library(rstan)
library(stringr)
library(ggh4x)
library(patchwork)

source("scripts/common/mcmc-util.R")
source("scripts/common/setup-argparse.R")
source("scripts/common/plot-settings.R")

parser$add_argument("--pf-data-list-format")
parser$add_argument("--pf-prior-optim-stan-model")
parser$add_argument("--output-pf-prior-plot")
args <- parser$parse_args()

list_data <- readRDS(args$pf_data_list_format)

n_icu_stays <- length(list_data)
n_prior_samples <- 50000

prefit <- stan_model(args$pf_prior_optim_stan_model)

res <- mclapply(1 : n_icu_stays, mc.cores = 5, function(icustay_index) {
  indiv_prior_data <- list_data[[icustay_index]]

  n_spline_coef <- indiv_prior_data %>%
    extract2('x_obs_mat') %>%
    ncol()

  ## prior samples occurr on the y-rescaled scale
  beta_zero_prior_samples <- rnorm(
    n = n_prior_samples,
    mean = 0,
    sd = 1
  )

  spline_coef_prior_samples <- array(
    data = NA,
    dim = c(n_prior_samples, n_spline_coef),
    dimnames = list(
      sprintf('iteration_%d', 1 : n_prior_samples),
      sprintf('zeta_%d', 1 : n_spline_coef)
    )
  )

  # first spline coef is indep
  spline_coef_prior_samples[, 1] <- rnorm(
    n = n_prior_samples,
    mean = 0,
    sd = 0.5
  )

  # rw smoothing prior on the rest
  for (coef_index in 2 : n_spline_coef) {
    spline_coef_prior_samples[, coef_index] <- rnorm(
      n = n_prior_samples,
      mean = spline_coef_prior_samples[, coef_index - 1],
      sd = 0.5
    )
  }

  event_samples <- array(
    data = NA,
    dim = c(n_prior_samples, 2),
    dimnames = list(
      sprintf('iteration_%d', 1 : n_prior_samples),
      c(
        sprintf('event_time[%d]', icustay_index),
        sprintf('event_indicator[%d]', icustay_index)
      )
    )
  )

  for (iteration_index in 1 : n_prior_samples) {
    beta_zero <- beta_zero_prior_samples[iteration_index]
    spline_coef <- spline_coef_prior_samples[iteration_index, ]

    boundary_knots <- indiv_prior_data$boundary_knots
    internal_knots <- indiv_prior_data$internal_knots
    scaled_threshold <- indiv_prior_data$threshold_ctr
    init <- boundary_knots[1] + 1e-4

    eval_spline <- function(x) {
      sp_mat <- bSpline(x, knots = internal_knots, Boundary.knots = boundary_knots)
      res <- beta_zero + sp_mat %*% spline_coef
      return(res)
    }

    root_f <- function(x) {
      eval_spline(x) - scaled_threshold
    }

    # need to do some error handling for no root splines
    res_rootsolve <- rootSolve::uniroot.all(
      f = root_f,
      interval = c(
        max(boundary_knots[1], 0),
        boundary_knots[2]
      )
    )

    if (length(res_rootsolve) == 0) {
      # if there is no root, then we don't observe the respiratory failure
      # event and they are instead assumed to be discharged, or expire, at
      # their final observation.
      event_time <- boundary_knots[2]
      event_indicator <- 0
    } else {
      event_time <- min(res_rootsolve)

      if (event_time == 0) {
        event_time = 1e-8 # otherwise the beta distribution falls over
      }

      # this is an annoying numeric edge case
      if (event_time == boundary_knots[2]) {
        event_indicator <- 0
      } else {
        event_indicator <- 1
      }
    }

    event_samples[iteration_index, ] <- c(event_time, event_indicator)
  }

  return(event_samples)
})

parameter_res <- lapply(1 : n_icu_stays, function(icustay_index) {
  indiv_prior_data <- list_data[[icustay_index]]
  indiv_prior_samples <- res[[icustay_index]]

  stan_data <- list(
    n_prior_samples = n_prior_samples,
    event_time = indiv_prior_samples[, 1],
    event_indicator = indiv_prior_samples[, 2],
    lower_limit = indiv_prior_data$boundary_knots[1],
    upper_limit = indiv_prior_data$boundary_knots[2]
  )

  optm_list <- lapply(1 : 10, function(x) {
    optimizing(
      prefit,
      data = stan_data
    )
  })

  optm_best <- lapply(optm_list, function(x) {
    x$value
  }) %>%
    which.max()

  optm_res <- optm_list[[optm_best]]

  pars <- as.list(optm_res$par) %>%
    c(id = icustay_index) %>%
    as_tibble()

  # this function will need to be reused later -- pull into a file then.
  plot_f <- function(x) {
    optm_res$par['weight'] * dbeta(
      (x - stan_data$lower_limit)/ (stan_data$upper_limit - stan_data$lower_limit),
      shape1 = optm_res$par['beta_alpha'],
      shape2 = optm_res$par['beta_beta'],
    ) *
      (1 / (stan_data$upper_limit - stan_data$lower_limit)) # jacobian

    # There is a difficulty in visualising this 'density', because the histogram
    # normalisiation doesn't work (the spike / atom cannot be bigger than
    # one by construction, but the exact normalisation is not clear
  }

  plot_tbl <- tibble(
    x = seq(
      from = max(stan_data$lower_limit, 0) + 0.01,
      to = stan_data$upper_limit - 0.01,
      length.out = 500
    ),
    id = icustay_index
  )  %>%
    mutate(y = plot_f(x))

  inner_res <- list(
    pars = pars,
    plot_tbl = plot_tbl
  )

})

param_tbl <- bind_named_sublists(parameter_res, 'pars', 1) %>%
  as_tibble() %>%
  mutate(
    plot_label = sprintf(
      'pi = %.5f, a = %.2f, b = %.2f',
      weight,
      beta_alpha,
      beta_beta
    )
  )

sample_tbl <- lapply(res, function(sub_list) {
  id <- str_extract(names(sub_list[1, ])[1], '(\\d+)') %>%
    as.integer()

  tibble(
    event_time = sub_list[, 1],
    event_indicator = as.integer(sub_list[, 2]),
    id = id
  )
}) %>%
  bind_rows() %>%
  mutate(
    plot_type = 'histogram'
  ) %>%
  left_join(param_tbl, by = 'id')

plot_tbl <- bind_named_sublists(parameter_res, 'plot_tbl', 1) %>%
  as_tibble() %>%
  left_join(param_tbl, by = 'id') %>%
  mutate(plot_type = 'fitted_dens')

plot_list <- lapply(1 : n_icu_stays, function(icustay_index) {
  sub_plot_tbl <- plot_tbl %>%
    filter(id == icustay_index)

  sub_sample_tbl <- sample_tbl %>%
    filter(id == icustay_index)

  bw <- sub_sample_tbl %>%
    filter(event_indicator == 1) %>%
    pull(event_time) %>%
    bw.SJ()

  weight_val <- sub_sample_tbl %>%
    pull(weight) %>%
    unique()

  ggplot(data = sub_plot_tbl) +
    geom_line(aes(x = x, y = y)) +
    geom_histogram(
      data = sub_sample_tbl %>% filter(event_indicator == 1),
      mapping = aes(x = event_time, y = after_stat(density * weight_val)),
      alpha = 0.5,
      binwidth = bw,
      fill = blues[2]
    ) +
    geom_col(
      data = sub_sample_tbl %>%
        filter(event_indicator == 0) %>%
        distinct(),
      mapping = aes(x = event_time, y = (1 - weight)),
      alpha = 0.7,
      width = bw,
      fill = highlight_col
    ) +
    xlab(bquote("T"[.(icustay_index)])) +
    ylab(bquote("p"[1,.(icustay_index)]("T"[.(icustay_index)]))) +
    ggtitle(label = bquote(italic('i')==.(icustay_index)))
})

p1 <- wrap_plots(plot_list, ncol = 3)

ggsave_base(
  filename = args$output_pf_prior_plot,
  plot = p1,
  height = 80,
  width = 20
)

p2 <- wrap_plots(plot_list %>% inset(1 : (n_icu_stays - 10), NULL), ncol = 3)

ggsave_fullpage(
  filename = str_replace(args$output_pf_prior_plot, '.png', '-small.png'),
  plot = p2,
)

saveRDS(
  file = args$output,
  object = param_tbl
)
