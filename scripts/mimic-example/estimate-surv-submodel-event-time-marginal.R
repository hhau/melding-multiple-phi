library(dplyr)
library(rstan)
library(magrittr)

source('scripts/common/logger-setup.R')
source('scripts/common/setup-argparse.R')
source("scripts/common/plot-settings.R")
source('scripts/common/mcmc-util.R')

parser$add_argument("--surv-prior-samples-raw")
parser$add_argument("--pf-data-list")
parser$add_argument("--pf-prior-optim-stan-model")
parser$add_argument("--output-surv-event-time-only-prior-plot")
args <- parser$parse_args()

raw_samples <- readRDS(args$surv_prior_samples_raw)
pf_list_data <- readRDS(args$pf_data_list)

flog.info(
  "MIMC-example-: estimating \\hat(p)_{2}(phi_12) parameters",
  name = base_filename
)

pf_prefit <- stan_model(args$pf_prior_optim_stan_model)
n_icu_stays <- length(raw_samples)

parameter_res <- lapply(1 : n_icu_stays, function(icustay_index) {
  indiv_prior_data <- pf_list_data[[icustay_index]]
  indiv_prior_samples <- raw_samples[[icustay_index]] %>%
    filter(complete.cases(.))

  stan_data <- list(
    n_prior_samples = nrow(indiv_prior_samples),
    event_time = indiv_prior_samples$event_time,
    event_indicator = indiv_prior_samples$event_indicator,
    lower_limit = 0,
    upper_limit = indiv_prior_data$boundary_knots[2]
  )

  optm_list <- lapply(1 : 10, function(x) {
    optimizing(
      pf_prefit,
      data = stan_data
    )
  })

  optm_best <- lapply(optm_list, function(x) {
    if (x$return_code != 0) {
      return(NA)
    } else {
      x$value
    }
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

sample_tbl <- lapply(raw_samples, function(sub_list) {
  local_tbl <- sub_list %>%
    filter(complete.cases(.))

  tibble(
    event_time = local_tbl$event_time,
    event_indicator = local_tbl$event_indicator,
    id = local_tbl$icustay_index
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
    ylab(bquote("p"[2 * ',' ~ .(icustay_index)]("T"[.(icustay_index)]))) +
    ggtitle(label = bquote(italic('i')==.(icustay_index)))
})

p1 <- wrap_plots(plot_list, ncol = 3)

ggsave_base(
  filename = args$output_surv_event_time_only_prior_plot,
  plot = p1,
  height = 80,
  width = 20
)

interesting_plot_ids <- c(6, 13, 15)

p2 <- wrap_plots(
  plot_list[interesting_plot_ids],
  ncol = 3
)

ggsave_fullpage(
  filename = str_replace(args$output_surv_event_time_only_prior_plot, '.png', '-small.pdf'),
  adjust_height = -15,
  plot = p2,
)

saveRDS(
  file = args$output,
  object = param_tbl
)
