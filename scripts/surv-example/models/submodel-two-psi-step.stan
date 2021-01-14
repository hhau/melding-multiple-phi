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
}

parameters {
  // intercept
  real beta_zero;

  // baseline coefficient
  real beta_one;

  // baseline hazard parameter(s) (weibull lambda)
  real <lower = 0> hazard_lambda;
}

model {
  // wiebull model, no censoring.
  // log harzard
  target += log(hazard_lambda);
  target += sum((hazard_lambda - 1) * log(event_time));
  target += beta_zero + sum(baseline_measurement * beta_one);

  // log survival probability
  target += - sum(
    pow_vec(event_time, hazard_lambda) .* 
    exp(beta_zero + baseline_measurement * beta_one)
  );

  // priors
  target += normal_lpdf(beta_zero | 0, 1);
  target += normal_lpdf(beta_one | 0, 1);
  target += normal_lpdf(hazard_lambda | 0, 1);
}
