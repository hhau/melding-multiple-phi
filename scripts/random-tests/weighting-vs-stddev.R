curve(
  dnorm,
  from = -5, 
  to = 5
)

q_w <- function(x) {
  dnorm(x) ^ 0.25
}

res <- integrate(q_w, lower = -15, upper = 15)

f_w <- function(x) {
  (dnorm(x) ^ 0.25) / res$value
}

curve(
  f_w,
  from = -5, 
  to = 5,
  col = "blue",
  add = T
)

wide_f <- function(x) {
  dnorm(x, sd = sqrt(1 / 0.25))
}

curve(
  wide_f,
  from = -5, 
  to = 5,
  col = "red",
  lty = "dashed",
  add = T
)


