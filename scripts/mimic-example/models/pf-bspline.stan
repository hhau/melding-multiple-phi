data {
  int <lower = 1> n_icu_stays;
  int <lower = 1> n_total_obs;
  int <lower = 1> n_plot_points;
  int <lower = 3> n_basis_coef;
  int <lower = 1, upper = n_total_obs + 1> subset_vector [n_icu_stays + 1];
  matrix [n_total_obs, n_basis_coef] obs_matrix;
  vector [n_total_obs] y_vec_ctr;
  vector [n_icu_stays] ctr_means;
  vector [n_icu_stays] ctr_sds;
  matrix [n_plot_points, n_basis_coef] plot_matrix [n_icu_stays];
  real <lower = 0> spline_coef_prior_sd;
  real <lower = 0> noise_df;
}

parameters {
  real <lower = 0> y_sigma;
  vector [n_icu_stays] beta_zero;
  vector [n_basis_coef] spline_coef [n_icu_stays];
}

model {
  for (ii in 1 : n_icu_stays) {
    // build each persons mean
    int n_obs_per_icu_stay = subset_vector[ii + 1] - subset_vector[ii];
    vector [n_obs_per_icu_stay] indiv_obs_mu;
    vector [n_obs_per_icu_stay] indiv_y_ctr;
    matrix [n_obs_per_icu_stay, n_basis_coef] indiv_obs_matrix;
    indiv_y_ctr = y_vec_ctr[subset_vector[ii] : (subset_vector[ii + 1] - 1)];

    indiv_obs_matrix = obs_matrix[
      subset_vector[ii] : (subset_vector[ii + 1] - 1),
    ];

    indiv_obs_mu = beta_zero[ii] + indiv_obs_matrix * spline_coef[ii];

    // add on the likelihood
    target += student_t_lpdf(indiv_y_ctr | noise_df, indiv_obs_mu, y_sigma);

    // spline coefficient prior for appropriate wigglyness.
    target += normal_lpdf(spline_coef[ii][1] | 0, spline_coef_prior_sd);

    for (jj in 2 : n_basis_coef) {
      real spline_delta = spline_coef[ii][jj] - spline_coef[ii][jj - 1];
      target += normal_lpdf(spline_delta | 0, spline_coef_prior_sd);
    }
  }

  // priors - response has been centred and scaled
  target += normal_lpdf(beta_zero | 0, 1);
  target += normal_lpdf(y_sigma | 0, 1);
}

generated quantities {
  matrix [n_icu_stays, n_plot_points] plot_mu;

  for (ii in 1 : n_icu_stays) {
    // the plot_mus are on the original scale.
    plot_mu[ii, ] = (ctr_means[ii] +
      (beta_zero[ii] + plot_matrix[ii] * spline_coef[ii]) * ctr_sds[ii])';
  }
}
