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
  // note that event_indicator will be a vector of ones n zeros, but 
  // we can't put the bounds on it because they will have log_prob of -inf
  vector <lower = 0> [n_patients] event_time;
  vector [n_patients] event_indicator;

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
  // Weibull model log hazard
  // note that this model is never used to compute a gradient
  // which is why I can do this branching 
  target += sum(event_indicator) * log(hazard_gamma);

  for (ii in 1 : n_patients) {
    if (event_indicator[ii] == 1) {
      target += (hazard_gamma - 1) * log(event_time[ii]);
      target += z_common[ii];
    }
  }

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
