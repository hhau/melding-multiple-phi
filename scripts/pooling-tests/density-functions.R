library(mvtnorm)
library(cubature)
suppressPackageStartupMessages(library(sm))

f_normal_log_generator <- function(outer_lambda_1, outer_lambda_2) {
  q_normal_log <- function(
    phi,
    lambda_1,
    lambda_2
  ) {
    res <- 
      lambda_1 * dnorm(t(phi)[, 1, drop = FALSE], mean = mean_p1, log = T) + 
      lambda_2 * dmvnorm(
        x = t(phi),
        sigma = sigma_mat,
        log = T
      ) +
      lambda_1 * dnorm(t(phi)[, 2, drop = FALSE], mean = mean_p3, log = T)
    return(exp(res))
  }

  nc_normal_log <- cubintegrate(
    q_normal_log,
    lower = c(-15, -15),
    upper = c(15, 15),
    method = "hcubature",
    lambda_1 = outer_lambda_1,
    lambda_2 = outer_lambda_2,
    nVec = 62500
  )

  f_normal_log <- function(phi) {
    q_normal_log(
      phi,
      lambda_1 = outer_lambda_1,
      lambda_2 = outer_lambda_2
    ) / nc_normal_log$integral
  }

  return(f_normal_log)  
}

f_normal_linear_generator <- function(outer_lambda_1, outer_lambda_2) {
  q_normal_linear <- function(
    phi,
    lambda_1 = outer_lambda_1,
    lambda_2 = outer_lambda_2
  ) {
    res <-
      (
        lambda_1 * dnorm(t(phi)[, 1, drop = FALSE], mean = mean_p1, sd = 1) +
        lambda_2 * dnorm(t(phi)[, 1, drop = FALSE], mean = 0, sd = 1)  
      ) *
      (
        lambda_2 * dnorm(t(phi)[, 2, drop = FALSE], mean = 0, sd = 1) +
        lambda_1 * dnorm(t(phi)[, 2, drop = FALSE], mean = mean_p3, sd = 1)
      ) 

    return(res)
  }

  nc_normal_linear <- cubintegrate(
    q_normal_linear,
    lower = c(-15, -15),
    upper = c(15, 15),
    method = "hcubature",
    lambda_1 = outer_lambda_1,
    lambda_2 = outer_lambda_2,
    nVec = 62500
  )

  f_normal_linear <- function(phi) {
    q_normal_linear(
      phi,
      lambda_1 = outer_lambda_1,
      lambda_2 = outer_lambda_2
    ) / nc_normal_linear$integral
  }

  return(f_normal_linear)  
}


f_student_t_log_generator <- function(outer_lambda_1, outer_lambda_2) {
  q_student_t_log <- function(
    phi,
    lambda_1 = outer_lambda_1,
    lambda_2 = outer_lambda_2
  ) {
    res <- 
      lambda_1 * dt(x = t(phi)[, 1, drop = FALSE] - mean_p1, df = t_df, log = T) + 
      lambda_2 * dmvt(
        x = t(phi),
        sigma = sigma_mat,
        df = t_df,
        log = T
      ) +
      lambda_1 * dt(x = t(phi)[, 2, drop = FALSE] - mean_p3, df = t_df, log = T)
    return(exp(res))
  }

  nc_student_t_log <- cubintegrate(
    q_student_t_log,
    lower = c(-15, -15),
    upper = c(15, 15),
    method = "hcubature",
    lambda_1 = outer_lambda_1,
    lambda_2 = outer_lambda_2,
    nVec = 62500
  )

  f_student_t_log <- function(phi) {
    q_student_t_log(
      phi,
      lambda_1 = outer_lambda_1,
      lambda_2 = outer_lambda_2
    ) / nc_student_t_log$integral
  }

  return(f_student_t_log)  
}

f_student_t_linear_generator <- function(outer_lambda_1, outer_lambda_2) {
  twod_samples <- rmvt(
    n = 2e4,
    sigma = sigma_mat,
    df = t_df
  )
  
  bandwidths <- apply(twod_samples, 2, bw.SJ) * 4

  kde_marginal <- function(phi, margin) {
    sm.density(
      x = twod_samples[, margin],
      h = bandwidths[margin],
      eval.points = as.vector(phi),
      display = "none"
    )[['estimate']]
  }

  q_student_t_linear <- function(
    phi,
    lambda_1 = outer_lambda_1,
    lambda_2 = outer_lambda_2
  ) {
    res <-
      (
        lambda_1 * dt(t(phi)[, 1, drop = FALSE] - mean_p1, df = t_df) +
        lambda_2 * kde_marginal(t(phi)[, 1, drop = FALSE], margin = 1)  
      ) *
      (
        lambda_2 * kde_marginal(t(phi)[, 2, drop = FALSE], margin = 2) +
        lambda_1 * dt(t(phi)[, 2, drop = FALSE] - mean_p3, df = t_df)
      ) 

  }
  
  nc_student_t_linear <- cubintegrate(
    q_student_t_linear,
    lower = c(-15, -15),
    upper = c(15, 15),
    method = "hcubature",
    lambda_1 = outer_lambda_1,
    lambda_2 = outer_lambda_2,
    nVec = 62500,
    maxEval = 1e3
  )

  f_student_t_linear <- function(phi) {
    q_student_t_linear(
      phi,
      lambda_1 = outer_lambda_1,
      lambda_2 = outer_lambda_2
    ) / nc_student_t_linear$integral
  }
  
}
