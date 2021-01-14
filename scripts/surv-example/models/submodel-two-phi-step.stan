data {
  int <lower = 1> n_patients;  
  vector [n_patients] baseline_measurement;

  // in the phi-step these are all fixed.
  vector [n_covariates] beta_zero;
  real mu_beta_zero;
  real <lower = 0> sigma_beta_zero;

  // baseline hazard parameter(s) (weibull lambda)
  vector <lower = 0> hazard_lambda;
}

parameters {
  // assuming no censoring for the moment
  // in phi-step this is a parameter
  real <lower = 0> event_time [n_patients];
}

model {
  // wiebull model, no censoring.
  target += log(hazard_lambda);
  target += (hazard_lambda - 1) * log(event_time);
  target += (baseline_measurement') * beta_zero;

  // priors
  target += normal_lpdf(beta_zero | mu_beta_zero, sigma_beta_zero);
  target += normal_lpdf(mu_beta_zero | 0, 1);
  target += normal_lpdf(sigma_beta_zero | 0, 1);
}