data {
  int <lower = 1> n_obs;
  int <lower = 1> n_patients;
  vector [n_obs] Y;
  int obs_ids [n_obs]; // because we lack ragged arrays
  vector [n_obs] obs_times; // times
  vector [n_obs] X; // covariate values

  // add prediction points? Maybe do in R?
  real y_threshold;

  // plotting x-values, this makes the samples array really big
  // but it is a lot easier to generate them here.
  int n_plot;
  vector [n_plot] t_plot;
  vector [n_plot] x_plot;
}

parameters {
  vector [n_patients] alpha_zero;
  vector [n_patients] alpha_one;
  vector [n_patients] beta_one;
  // real mu_alpha_zero;
  // real mu_alpha_one;
  // real mu_beta_one;
  // real <lower = 0> sigma_alpha_zero;
  // real <lower = 0> sigma_alpha_one;
  // real <lower = 0> sigma_beta_one;

  real <lower = 0> sigma_noise_x;
  real <lower = 0> sigma_noise_t;
}

transformed parameters {
  real <lower = 0> sigma_total = sqrt(sigma_noise_x^2 + sigma_noise_t ^2);
}

model {
  // dgp
  vector [n_obs] temp_mu;

  for (ii in 1 : n_obs) {
    temp_mu[ii] = 
      alpha_zero[obs_ids[ii]] + 
      alpha_one[obs_ids[ii]] * X[ii] +
      beta_one[obs_ids[ii]] * obs_times[ii];
  }

  target += normal_lpdf(Y | temp_mu, sigma_total);

  // priors
  // target += normal_lpdf(alpha_zero | mu_alpha_zero, sigma_alpha_zero);
  // target += normal_lpdf(alpha_one | mu_alpha_one, sigma_alpha_one);
  // target += normal_lpdf(beta_one | mu_beta_one, sigma_beta_one);

  target += normal_lpdf(alpha_zero | 1.0, 0.1);
  target += normal_lpdf(alpha_one | 0.5, 0.1);
  target += log_mix(
    0.5,
    normal_lpdf(beta_one | -1.0, 0.2),
    normal_lpdf(beta_one | 0.0, 0.1)
  );

  // target += normal_lpdf(mu_alpha_zero | 1.0, 1.0);
  // target += normal_lpdf(mu_alpha_one | 0.5, 1.0);
  // target += normal_lpdf(mu_beta_one | 0.0, 1.0);

  // target += normal_lpdf(sigma_alpha_zero | 1.0, 0.5);
  // target += normal_lpdf(sigma_alpha_one | 1.0, 0.5);
  // target += normal_lpdf(sigma_beta_one | 1.0, 0.5);
  
  target += normal_lpdf(sigma_noise_x | 0.2, 0.2);
  target += normal_lpdf(sigma_noise_t | 0.02, 0.02);
}

generated quantities {
  // Mon  3 Aug 10:58:03 2020 - not sure this is right now that I have
  // the regression aspect of the model? 
  vector [n_patients] event_time = (y_threshold - alpha_zero) ./ beta_one;
  matrix [n_patients, n_plot] plot_mu = rep_matrix(alpha_zero, n_plot) + alpha_one * (x_plot') + beta_one * (t_plot'); 
}
