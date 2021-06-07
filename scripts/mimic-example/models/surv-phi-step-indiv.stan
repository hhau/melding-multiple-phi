data {
  int <lower = 1> n_theta; // includes intercept
  int <lower = 1> n_segments;
}

parameters {
  // need to include the baseline measurements and thetas because the are
  // multiplicative modifiers of the survival probability and thus do not cancel
  // in the log-likelihood.
  vector [n_theta] theta;
  row_vector [n_theta] baseline_data_x;

  real breakpoint_lower;
  real breakpoint_upper;

  // baseline hazard parameter(s) (Weibull gamma)
  real <lower = 0> hazard_gamma;

  // longitudinal associative strength (alpha)
  real alpha;

  // phi_{1 \cap 2} -- from the event time submodel
  real event_indicator;
  real <lower = 0> event_time;

  // phi_{2 \cap 3} -- from the longitudinal submodel
  // used to calculate m'_{i}(t)
  real breakpoint;
  vector [n_segments] eta_slope; // first segment is eta^{b} -- eta 'before'
}

transformed parameters {
  real width = breakpoint_upper - breakpoint_lower;
  real breakpoint_raw = (breakpoint - breakpoint_lower) / width;
}

model {
  real temp_surv_prob_common;

  if (event_indicator == 1) { // add the hazard
    target += log(hazard_gamma) + (hazard_gamma - 1) * log(event_time);
    target += baseline_data_x * theta; // constant, but easier to leave in

    if (event_time < breakpoint) {
      target += alpha * eta_slope[1];
    } else {
      target += alpha * eta_slope[2];
    }
  }

  // always add the survival probability
  // this term is common despite to both cases
  temp_surv_prob_common = -exp(baseline_data_x * theta);

  if (event_time < breakpoint) {
    real t1 = exp(alpha * eta_slope[1]);
    real t2 = pow(event_time, hazard_gamma);
    real temp_lower =  t1 * t2;
    target += temp_surv_prob_common * temp_lower;
  } else {
    real t1 = exp(alpha * eta_slope[1]);
    real t2 = pow(breakpoint, hazard_gamma);
    real t3 = exp(alpha * eta_slope[2]);
    real t4 = pow(event_time, hazard_gamma);
    real t5 = pow(breakpoint, hazard_gamma);
    real temp_upper = (t1 * t1) + (t3) * (t4 - t5);
    target += temp_surv_prob_common * temp_upper;
  }

  // priors for phi_{1 \cap 2} (implicit) and phi_{2 \cap 3} (explicit)
  // align priors with the fluid submodel
  target += gamma_lpdf(eta_slope | 1.53, 0.24);

  // prior for kappa + jacobian for linear change of variables
  target += beta_lpdf(breakpoint_raw | 5.0, 5.0);
  target += -log(width);
}
