library(mvtnorm)
library(cubature)

f_normal_phi_12_marginal <- function(phi_12, log_scale = FALSE) {
  dnorm(phi_12, mean = mean_p1, sd = 1, log = log_scale)
}

f_normal_phi_23_marginal <- function(phi_23, log_scale = FALSE) {
  dnorm(phi_23, mean = mean_p3, sd = 1, log = log_scale)
}

f_normal_twod_marginal <- function(phi, log_scale = FALSE) {
  dmvnorm(
    x = phi,
    mean = rep(0, 2),
    sigma = sigma_mat
  )
} 

f_normal_logarithmic_generator <- function(outer_lambda_1, outer_lambda_2) {
  q_normal_log <- function(
    phi,
    lambda_1,
    lambda_2
  ) {
    res <- 
      lambda_1 * dnorm(t(phi)[, 1, drop = FALSE], mean = mean_p1, log = TRUE) + 
      lambda_2 * dmvnorm(
        x = t(phi),
        sigma = sigma_mat,
        log = TRUE
      ) +
      lambda_1 * dnorm(t(phi)[, 2, drop = FALSE], mean = mean_p3, log = TRUE)
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

