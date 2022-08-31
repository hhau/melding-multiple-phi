#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix spline_basis_inside_rcpp(NumericVector x) {
  Environment splines2_pkg = Environment::namespace_env("splines2");
  Function ispline_basis = splines2_pkg["iSpline"];
  NumericVector knots_vec = {0.5};
  NumericVector boundary_knots = {0.0, 1.0};
  
  NumericMatrix result = ispline_basis(
    x,
    Named("knots", knots_vec),
    Named("Boundary.knots", boundary_knots)
  );
  
  return result;
}
