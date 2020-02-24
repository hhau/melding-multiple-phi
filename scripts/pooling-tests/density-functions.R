library(mvtnorm)
library(cubature)

f_normal_log_generator <- function(outer_w_1, outer_w_2, outer_w_3) {
  
  q_normal_log <- function(
    phi,
    w_1 = outer_w_1,
    w_2 = outer_w_2,
    w_3 = outer_w_3
  ) {
    res <- 
      w_1 * dnorm(phi[1], mean = -2, sd = 1, log = T) + 
      w_2 * dmvnorm(
        x = as.vector(phi),
        sigma = matrix(
          c(1, 0.8, 0.8, 1),
          nrow = 2,
          ncol = 2,
          byrow = TRUE
        ),
        log = T
      ) +
      w_3 * dnorm(phi[2], mean = 2, sd = 1, log = T)
    return(exp(res))
  }

  nc_normal_log <- cubintegrate(
    q_normal_log,
    lower = c(-15, -15),
    upper = c(15, 15),
    method = "hcubature",
    w_1 = outer_w_1,
    w_2 = outer_w_2,
    w_3 = outer_w_3
  )

  f_normal_log <- function(phi) {
    q_normal_log(
      phi,
      w_1 = outer_w_1,
      w_2 = outer_w_2,
      w_3 = outer_w_3
    ) / nc_normal_log$integral
  }

  return(f_normal_log)  
}

f_student_t_log_generator <- function(outer_w_1, outer_w_2, outer_w_3) {
  
  q_student_t_log <- function(phi,
    w_1 = outer_w_1,
    w_2 = outer_w_2,
    w_3 = outer_w_3
  ) {
    res <- 
      w_1 * dt(x = phi[1] + 2, df = 4, log = T) + 
      w_2 * dmvt(
        x = as.vector(phi),
        sigma = matrix(
          c(1, 0.8, 0.8, 1),
          nrow = 2,
          ncol = 2,
          byrow = TRUE
        ),
        df = 4,
        log = T
      ) +
      w_3 * dt(x = phi[2] - 2, df = 4, log = T)
    return(exp(res))
  }

  nc_student_t_log <- cubintegrate(
    q_student_t_log,
    lower = c(-15, -15),
    upper = c(15, 15),
    method = "hcubature",
    w_1 = outer_w_1,
    w_2 = outer_w_2,
    w_3 = outer_w_3
  )

  f_student_t_log <- function(phi) {
    q_student_t_log(
      phi,
      w_1 = outer_w_1,
      w_2 = outer_w_2,
      w_3 = outer_w_3
    ) / nc_student_t_log$integral
  }

  return(f_student_t_log)  
}

