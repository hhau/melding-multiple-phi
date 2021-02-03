functions {
  vector pow_vec(vector x, real y) {
    int N = rows(x);
    vector [N] res;
    for (nn in 1 : N) {
      res[nn] = pow(x[nn], y);
    }
    return res;
  }
}

data {
  // global constant
  int <lower = 1> n_patients;
  int <lower = 1> n_long_beta;

  // Y_{2} - baseline data
  vector [n_patients] baseline_measurement;
  real log_crude_event_rate; // technically comes from Y_{1} but eh.
}

parameters {
  // intercept
  real beta_zero;

  // baseline coefficient
  real beta_one;

  // baseline hazard parameter(s) (Weibull gamma)
  real <lower = 0> hazard_gamma;

  // longitudinal associative strength (alpha)
  real alpha;

  // phi_{1 \cap 2} -- from the event time submodel
  // assuming no censoring for the moment
  vector <lower = 0> [n_patients] event_time;

  // phi_{2 \cap 3} -- from the longitudinal submodel
  // used to calculate m_{i}(t)
  matrix [n_patients, n_long_beta] long_beta;
}

transformed parameters {
  vector [n_patients] z_common;

  z_common = beta_zero + 
    baseline_measurement * beta_one + 
    alpha * long_beta[, 1];
}

model {
  // Weibull model, no censoring.
  // log hazard
  target += n_patients * log(hazard_gamma);
  target += (hazard_gamma - 1) * sum(log(event_time));
  target += sum(z_common);

  // log survival probability
  target += -sum(
    pow_vec(event_time, hazard_gamma) .* z_common
  );

  // priors 
  target += normal_lpdf(beta_zero | log_crude_event_rate, 1);
  target += normal_lpdf(beta_one | 0, 1);
  target += normal_lpdf(hazard_gamma | 0, 1);
  target += normal_lpdf(alpha | 0, 1);
}