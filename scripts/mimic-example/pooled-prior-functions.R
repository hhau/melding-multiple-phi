library(dplyr)
library(mvtnorm)

pf_list_data <- readRDS('rds/mimic-example/submodel-1-pf-data-list-format.rds')
pf_prior_params <- readRDS('rds/mimic-example/submodel-1-marginal-prior-parameter-estimates.rds')
surv_prior_params <- readRDS('rds/mimic-example/submodel-2-marginal-prior-parameter-estimates.rds')
fluid_stan_data <- readRDS('rds/mimic-example/submodel-3-fluid-data-stan-format.rds')

# we know we have to do individual-at-a-time updating, so lets write the
# functions with that in mind.

# the ordering is (phi_12) = (event_time, event_indicator)
# (phi_23) = (breakpoint, eta_before, eta_after) -- as in the stan model.
# I need to check that this aligns with the second stage

submodel_1_prior_marginal <- function(phi_12, id, return_log = TRUE) {
  indiv_beta_pars <- pf_prior_params[id, ]
  indiv_data <- pf_list_data[[id]]
  lower_limit <- indiv_data$boundary_knots[1]
  upper_limit <- indiv_data$boundary_knots[2]
  length_of_stay <- upper_limit - lower_limit

  # flog.info(
  #   sprintf(
  #     "p_{1}(phi_12), id: %d, phi_12[1] = %f, phi_12[2] = %d, lower_limit: %f, upper_limit: %f",
  #     id, phi_12[1], phi_12[2], lower_limit, upper_limit
  #   ),
  #   name = base_filename
  # )

  stopifnot(
    (phi_12[1] >= lower_limit),
    (phi_12[1] <= upper_limit),
    (phi_12[2] == 1) | (phi_12[2] == 0)
  )

  scaled_event_time <- (phi_12[1] - lower_limit) / length_of_stay

  if (scaled_event_time == 0) {
    scaled_event_time <- 1e-8
  }

  if ((scaled_event_time == 1) & (phi_12[2] == 1)) {
    scaled_event_time <- 1 - 1e-8
  }

  with(indiv_beta_pars, {
    if (phi_12[2] == 1) {
      lp <- log(weight) + dbeta(scaled_event_time, beta_alpha, beta_beta, log = TRUE)
    } else {
      lp <- log(1 - weight)
    }

    lp <- lp - log(length_of_stay)

    ifelse(return_log, return(lp), return(exp(lp)))
  })
}

logit_log_jac <- function(x, lb, ub) {
  log((1 / (ub - x)) + (1 / (x - lb)))
}

submodel_2_prior_marginal <- function(phi_12, phi_23, id, return_log = TRUE) {
  indiv_params_data <- surv_prior_params[[id]]
  with(indiv_params_data, {
    # flog.info(
    #   sprintf(
    #     "p_{2}(phi_12, phi_23), id: %d,
    #     phi_12[1] = %f, phi_12[2] = %d, lower_event_time_limit = %f, upper_event_time_limit = %f
    #     phi_23[1] = %f, phi_23[2] = %f, phi_23[3] = %f, lower_breakpoint_limit = %f, upper_breakpoint_limit = %f",
    #     id,
    #     phi_12[1], phi_12[2], lower_event_time_limit, upper_event_time_limit,
    #     phi_23[1], phi_23[2], phi_23[3], lower_breakpoint_limit, upper_breakpoint_limit
    #   ),
    #   name = base_filename
    # )

    stopifnot(
      (phi_12[1] >= lower_event_time_limit),
      (phi_12[1] <= upper_event_time_limit),
      (phi_12[2] == 1) | (phi_12[2] == 0),
      (phi_23[1] > lower_breakpoint_limit),
      (phi_23[1] < upper_breakpoint_limit),
      (phi_23[2] > 0),
      (phi_23[3] > 0)
    )

    if (phi_12[1] == lower_event_time_limit) {
      phi_12[1] <- lower_event_time_limit + 1e-8
    }

    if (phi_12[1] == upper_event_time_limit & phi_12[2] == 1) {
      phi_12[1] <- upper_event_time_limit - 1e-8
    }

    unconstrained_event_time <- suppressWarnings(brms::logit_scaled(
      x = phi_12[1],
      lb = lower_event_time_limit,
      ub = upper_event_time_limit
    ))

    unconstrained_breakpoint <- brms::logit_scaled(
      x = phi_23[1],
      lb = lower_breakpoint_limit,
      ub = upper_breakpoint_limit
    )

    unconstrained_eta_before <- log(phi_23[2])
    unconstrained_eta_after <- log(phi_23[3])

    if (phi_12[2] == 1) {
      par_vec <- c(
        unconstrained_event_time,
        unconstrained_breakpoint,
        unconstrained_eta_before,
        unconstrained_eta_after
      )

      lp <- dmvnorm(
        x = par_vec,
        mean = indiv_params_data$rf_event_mu,
        sigma = indiv_params_data$rf_event_sigma_mat,
        log = TRUE,
        checkSymmetry = FALSE
      )

      lp <- lp + unname(log(mix_weight))

      log_jac <- logit_log_jac(
        phi_12[1],
        lower_event_time_limit,
        upper_event_time_limit
      ) +
        logit_log_jac(
          phi_23[1],
          lower_breakpoint_limit,
          upper_breakpoint_limit
        ) -
        log(phi_23[2]) -
        log(phi_23[3])

    } else {
      par_vec <- c(
        unconstrained_breakpoint,
        unconstrained_eta_before,
        unconstrained_eta_after
      )

      lp <- dmvnorm(
        x = par_vec,
        mean = indiv_params_data$dd_event_mu,
        sigma = indiv_params_data$dd_event_sigma_mat,
        log = TRUE,
        checkSymmetry = FALSE
      )

      lp <- lp + unname(log(1 - mix_weight))
      log_jac <- logit_log_jac(
        phi_23[1],
        lower_breakpoint_limit,
        upper_breakpoint_limit
      ) -
        log(phi_23[2]) -
        log(phi_23[3])
    }

    overall <- lp + log_jac
    ifelse(return_log, return(overall), return(exp(overall)))
  })
}

submodel_3_prior_marginal <- function(phi_23, id, return_log = TRUE) {
  breakpoint_lower <- fluid_stan_data$breakpoint_lower[id]
  breakpoint_upper <- fluid_stan_data$breakpoint_upper[id]
  breakpoint_width <- breakpoint_upper - breakpoint_lower

  stopifnot(
    phi_23[1] > breakpoint_lower,
    phi_23[1] < breakpoint_upper,
    phi_23[2] > 0,
    phi_23[3] > 0
  )

  breakpoint_raw <- (phi_23[1] - breakpoint_lower) / breakpoint_width

  lp <- dbeta(breakpoint_raw, 5.0, 5.0, log = TRUE) +
    sum(dgamma(phi_23[2 : 3], shape = 1.53, rate = 0.24, log = TRUE))

  log_jac <- -log(breakpoint_width)
  overall <- lp + log_jac
  ifelse(return_log, return(overall), return(exp(overall)))
}

linear_pooled_prior <- function(phi_12, phi_23, lambda, id, return_log = TRUE) {
  # TODO: requires estimating p_{2}(phi_12) and p_{2}(phi_23), which is an effort
  # I am not ready to make now
  # Also eliminates correlation, so not particularly interesting tbh.
  stop('linear_pooled_prior not yet implemented')
}

logarithmic_pooled_prior <- function(phi_12, phi_23, lambda, id, return_log = TRUE) {
  stopifnot(length(lambda) == 3)
  lp1 <- lambda[1] * submodel_1_prior_marginal(phi_12, id)
  lp2 <- lambda[2] * submodel_2_prior_marginal(phi_12, phi_23, id)
  lp3 <- lambda[3] * submodel_3_prior_marginal(phi_23, id)

  overall <- lp1 + lp2 + lp3
  ifelse(return_log, return(overall), return(exp(overall)))
}

logarithmic_pooled_prior_overall <- function(phi_12, phi_23, lambda, id, return_log = TRUE) {
  lp <- logarithmic_pooled_prior(phi_12, phi_23, lambda, id) -
    submodel_1_prior_marginal(phi_12, id) -
    submodel_2_prior_marginal(phi_12, phi_23, id) -
    submodel_3_prior_marginal(phi_23, id)

  ifelse(return_log, return(lp), return(exp(lp)))
}
