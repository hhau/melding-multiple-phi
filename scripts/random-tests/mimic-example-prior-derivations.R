library(rriskDistributions)

# for the y-intercept in the cumulative fluid model
get.lnorm.par(
  p = c(0.025, 0.5, 0.99),
  q = c(500 / 1000, 5000 / 1000, 15000 / 1000),
  tol = 0.001
)

hist(rlnorm(n = 1e6, meanlog = 1.61, sdlog = 0.47), breaks = 100)

# eh, bit light on the left hand tail for my liking, but close enough.
qlnorm(c(0.025, 0.5, 0.99), 1.61, 0.47)

# -----
# For the RF hazard shape parameter in the survival models
get.gamma.par(
  p = c(0.01, 0.5, 0.99),
  q = c(1 / 5, 1, 2)
)

hist(rgamma(n = 5e5, shape = 9.05, rate = 8.72), breaks = 500)

# slopes -- definitely can't be zero, but also can't be huge.
get.gamma.par(
  p = c(0.025, 0.5, 0.975),
  q = c(0.5, 5, 20)
)

hist(rgamma(n = 5e5, shape = 1.53, rate = 0.24), breaks = 500, freq = FALSE)
