library(splines2)
library(Rcpp)

n_internal_knot <- 1

r_f <- function(x) {
  iSpline(
    x = x, 
    knots = seq(
      from = 0, 
      to = 1, 
      length.out = n_internal_knot + 2
    )[-c(1, n_internal_knot + 2)],
    Boundary.knots = c(0, 1)
  )
}

sourceCpp("scripts/random-tests/call-ispline-within-Rcpp.cpp")

spline_basis_inside_rcpp(0.2)

bench::mark(
  r_ver = r_f(0.2),
  rcpp_ver = spline_basis_inside_rcpp(0.2)
)
