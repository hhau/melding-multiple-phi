data {
  int <lower = 1> n_obs;
  int <lower = 1> n_patients;
  vector [n_obs] Y;
  int obs_ids [n_obs]; // because we lack ragged arrays
  vector [n_obs] obs_times; // X

  int n_plot;
  vector [n_plot] x_plot;

  int n_long_beta;
}

parameters {
  vector [n_long_beta] mu_beta;
  vector <lower = 0> [n_long_beta] sigma_beta;
  matrix [n_patients, n_long_beta] beta;
  real <lower = 0> sigma_y;
}

model {
  // dgp
  vector [n_obs] temp_mu;

  for (ii in 1 : n_obs) {
    temp_mu[ii] = beta[obs_ids[ii], 1] + 
      beta[obs_ids[ii], 2] * obs_times[ii];
  }

  target += normal_lpdf(Y | temp_mu, sigma_y);

  // priors
  for (ii in 1 : n_long_beta) {
    target += normal_lpdf(beta[, ii] | mu_beta[ii], sigma_beta[ii]);
  }
  target += normal_lpdf(mu_beta | 1.0, 1.0);
  target += lognormal_lpdf(sigma_beta | 0.0, 1.0);
  target += lognormal_lpdf(sigma_y | 0.0, 1.0);
}

generated quantities {
  matrix [n_patients, n_plot] plot_mu = rep_matrix(beta[, 1], n_plot) + beta[, 2] * x_plot';
}
