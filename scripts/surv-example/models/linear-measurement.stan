data {
  int <lower = 1> n_obs;
  int <lower = 1> n_patients;
  vector [n_obs] Y;
  int obs_ids [n_obs]; // because we lack ragged arrays
  vector [n_obs] obs_times; // X

  // add prediction points? Maybe do in R?
  real y_threshold;

  int n_plot;
  vector [n_plot] x_plot;
}

parameters {
  vector [n_patients] beta_zero;
  vector [n_patients] beta_one;
  real mu_beta_zero;
  real mu_beta_one;
  real <lower = 0> sigma_beta_zero;
  real <lower = 0> sigma_beta_one;
  real <lower = 0> sigma_y;
}

model {
  // dgp
  vector [n_obs] temp_mu;

  for (ii in 1 : n_obs) {
    temp_mu[ii] = beta_zero[obs_ids[ii]] + beta_one[obs_ids[ii]] * obs_times[ii];
  }

  target += normal_lpdf(Y | temp_mu, sigma_y);

  // priors
  target += normal_lpdf(beta_zero | mu_beta_zero, sigma_beta_zero);
  target += normal_lpdf(beta_one | mu_beta_one, sigma_beta_one);
  target += normal_lpdf(mu_beta_zero | 1.0, 1.0);
  target += normal_lpdf(mu_beta_one | 0.0, 1.0);
  target += lognormal_lpdf(sigma_beta_zero | 0.0, 1.0);
  target += lognormal_lpdf(sigma_beta_one | 0.0, 1.0);
  target += lognormal_lpdf(sigma_y | 0.0, 1.0);
}

generated quantities {
  vector [n_patients] event_time = (y_threshold - beta_zero) ./ beta_one;
  matrix [n_patients, n_plot] plot_mu = rep_matrix(beta_zero, n_plot) + beta_one * (x_plot'); 
}
