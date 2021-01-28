functions {
  real integrand (real x, real xc, real[] theta, real[] x_r, int[] x_i) {
    real hazard_gamma = theta[1];
    real alpha = theta[2];
    real beta_one = theta[3];
    real res = x ^(hazard_gamma - 1) * exp(alpha * beta_one * x);
    return(res);
  }
}

data {
  real hazard_gamma;
  real alpha;
  real beta_one;
  real event_time;
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  real x;
}

model {
  x ~ std_normal();
}

generated quantities {
  real integral_res = integrate_1d(
    integrand,
    0.0,
    event_time,
    {hazard_gamma, alpha, beta_one},
    x_r,
    x_i,
    0.001
  );
}
