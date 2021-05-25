library(rriskDistributions)

# for the y-intercept in the cumulative fluid model
get.lnorm.par(
  p = c(0.025, 0.5, 0.99),
  q = c(500, 5000, 15000),
  tol = 0.001
)

hist(rlnorm(n = 1e6, meanlog = 8.52, sdlog = 0.47), breaks = 100)

# eh, bit light on the left hand tail for my liking, but close enough.
qlnorm(c(0.025, 0.5, 0.99), 8.52, 0.47)

# -----
# For the RF hazard shape parameter in the survival models
get.gamma.par(
  p = c(0.01, 0.5, 0.99),
  q = c(1 / 5, 1, 2)
)

hist(rgamma(n = 5e5, shape = 9.05, rate = 8.72), breaks = 500)
