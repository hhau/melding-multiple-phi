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
  int <lower = 1> n_patients;
  vector [n_patients] baseline_measurement;

  // assuming no censoring for the moment
  vector <lower = 0> [n_patients] event_time;
  real log_crude_event_rate;
}

parameters {
  // intercept
  real beta_zero;

  // baseline coefficient
  real beta_one;

  // baseline hazard parameter(s) (Weibull gamma)
  real <lower = 0> hazard_gamma;
}

transformed parameters {
  vector [n_patients] eta = beta_zero + baseline_measurement * beta_one;
  real sum_eta = sum(eta);
}

model {
  // Weibull model, no censoring.
  // log hazard
  target += n_patients * log(hazard_gamma);
  target += sum((hazard_gamma - 1) * log(event_time));
  target += sum_eta;

  // log survival probability
  target += - sum(
    pow_vec(event_time, hazard_gamma) .* 
    exp(beta_zero + baseline_measurement * beta_one)
  );

  // priors - with ~600 uncensored events, these shouldn't matter too much
  target += normal_lpdf(beta_zero | log_crude_event_rate, 1);
  target += normal_lpdf(beta_one | 0, 1);
  target += normal_lpdf(hazard_gamma | 0, 1);
}
