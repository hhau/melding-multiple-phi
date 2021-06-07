data {
  int <lower = 1> n_icu_stays;
  int <lower = 1> n_total_obs;
  int <lower = 1> n_plot_points;
  int <lower = 1, upper = n_total_obs + 1> subset_vector [n_icu_stays + 1];
  vector [n_total_obs] y_vec;
  vector [n_total_obs] x_vec;
  vector [n_icu_stays] breakpoint_lower;
  vector [n_icu_stays] breakpoint_upper;
  matrix [n_icu_stays, n_plot_points] x_mat_plot;
}

transformed data {
  vector [n_icu_stays] widths = breakpoint_upper - breakpoint_lower;
}

parameters {
  vector <lower = 0, upper = 1> [n_icu_stays] breakpoint_raw;
  vector <lower = 0> [2] eta_slope [n_icu_stays];
  vector <lower = 0> [n_icu_stays] eta_zero_raw;
  real <lower = 0> y_sigma;
}

transformed parameters {
  vector [n_icu_stays] breakpoint = breakpoint_raw .* widths + breakpoint_lower;
  vector [n_total_obs] mu;
  vector [n_icu_stays] eta_zero;

  for (ii in 1 : n_icu_stays) {
    int obs_lower = subset_vector[ii];
    int obs_upper = subset_vector[ii + 1];
    int n_obs_per_icu_stay = obs_upper - obs_lower;
    vector [n_obs_per_icu_stay] indiv_obs_mu;
    vector [n_obs_per_icu_stay] indiv_obs_x = x_vec[obs_lower : (obs_upper - 1)];
    vector [n_obs_per_icu_stay] indiv_obs_y = y_vec[obs_lower : (obs_upper - 1)];

    for (jj in 1 : n_obs_per_icu_stay) {
      if (indiv_obs_x[jj] < breakpoint[ii]) {
        indiv_obs_mu[jj] = eta_zero_raw[ii] + eta_slope[ii][1] * breakpoint[ii] + eta_slope[ii][1] * (indiv_obs_x[jj] - breakpoint[ii]);
      } else {
        indiv_obs_mu[jj] = eta_zero_raw[ii] + eta_slope[ii][1] * breakpoint[ii] + eta_slope[ii][2] * (indiv_obs_x[jj] - breakpoint[ii]);
      }
    }

    mu[obs_lower : (obs_upper - 1)] = indiv_obs_mu;
    eta_zero[ii] = eta_zero_raw[ii] +  eta_slope[ii][1] * breakpoint[ii];
  }
}

model {
  target += normal_lpdf(y_vec | mu, y_sigma);
  target += lognormal_lpdf(eta_zero_raw | 1.61, 0.47);
  target += beta_lpdf(breakpoint_raw | 5.0, 5.0);
  target += normal_lpdf(y_sigma | 0.0, 5.0);

  for (ii in 1 : n_icu_stays) {
    target += gamma_lpdf(eta_slope[ii] | 1.53, 0.24);
  }
}

generated quantities {
  matrix [n_icu_stays, n_plot_points] plot_mu;

  for (ii in 1 : n_icu_stays) {
    for (jj in 1 : n_plot_points) {
      real temp_mu;

      if (x_mat_plot[ii, jj] < breakpoint[ii]) {
        temp_mu = eta_zero[ii] + eta_slope[ii][1] * (x_mat_plot[ii, jj] - breakpoint[ii]);
      } else {
        temp_mu = eta_zero[ii] + eta_slope[ii][2] * (x_mat_plot[ii, jj] - breakpoint[ii]);
      }

      plot_mu[ii, jj] = temp_mu;
    }
  }
}
