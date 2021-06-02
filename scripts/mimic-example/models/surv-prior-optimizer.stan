data {
  int <lower = 1> n_prior_samples;
  real lower_event_time_limit;
  real upper_event_time_limit;
  real lower_breakpoint_limit;
  real upper_breakpoint_limit;
  vector <lower = lower_event_time_limit, upper = upper_event_time_limit> [n_prior_samples] event_time_samples;
  int <lower = 1, upper = 2> event_indicator_samples [n_prior_samples];
  vector <lower = lower_breakpoint_limit, upper = upper_breakpoint_limit> [n_prior_samples] breakpoint_samples;
  vector <lower = 0> [n_prior_samples] eta_before_samples;
  vector <lower = 0> [n_prior_samples] eta_after_samples;
}

transformed data {
  // this might return +inf for z = 1?
  vector [n_prior_samples] unconstrained_event_time_samples;
  vector [n_prior_samples] unconstrained_breakpoint_samples = logit(
    (breakpoint_samples - lower_breakpoint_limit) / (upper_breakpoint_limit - lower_breakpoint_limit)
  );
  vector [n_prior_samples] unconstrained_eta_before_samples = log(eta_before_samples);
  vector [n_prior_samples] unconstrained_eta_after_samples = log(eta_after_samples);

  // numerical hack needed to get stan to not throw an initialisation error
  for (ii in 1 : n_prior_samples) {
    if (event_time_samples[ii] == upper_event_time_limit) {
      unconstrained_event_time_samples[ii] = logit(
        (event_time_samples[ii] - 1e-8 - lower_event_time_limit) /
          (upper_event_time_limit - lower_event_time_limit)
      );
    } else {
      unconstrained_event_time_samples[ii] = logit(
        (event_time_samples[ii] - lower_event_time_limit) /
          (upper_event_time_limit - lower_event_time_limit)
      );
    }
  }
}

parameters {
  vector [4] rf_event_mu;
  vector [3] dd_event_mu;
  cov_matrix [4] rf_event_sigma_mat;
  cov_matrix [3] dd_event_sigma_mat;
  real <lower = 0, upper = 1> mix_weight;
}

model {
  // likelihood
  for (ii in 1 : n_prior_samples) {
    if (event_indicator_samples[ii] == 1) { // rf event
      vector [4] temp_y = to_vector({
        unconstrained_event_time_samples[ii],
        unconstrained_breakpoint_samples[ii],
        unconstrained_eta_before_samples[ii],
        unconstrained_eta_after_samples[ii]
      });

      target += multi_normal_lpdf(temp_y | rf_event_mu, rf_event_sigma_mat);
      target += log(mix_weight);
    } else if (event_indicator_samples[ii] == 2) { // dd event
      vector [3] temp_y = to_vector({
        unconstrained_breakpoint_samples[ii],
        unconstrained_eta_before_samples[ii],
        unconstrained_eta_after_samples[ii]
      });

      target += multi_normal_lpdf(temp_y | dd_event_mu, dd_event_sigma_mat);
      target += log(1 - mix_weight);
    }
  }

  // jacobians? Really long, don't write twice, so omit them here?
  // are they needed? I think maybe they are? No, they are all constant w.r.t
  // the things in the 'parameters' block (samples r data), but the jacobian
  // will be necessary when we come to evaluate the density.

}
