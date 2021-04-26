data {
  int <lower = 1> n_theta; // includes intercept
  int <lower = 1> n_segments;
}

parameters {
  // baseline hazard parameter(s) (Weibull gamma)
  real <lower = 0> hazard_gamma;

  // longitudinal associative strength (alpha)
  real alpha;

  // phi_{1 \cap 2} -- from the event time submodel
  real event_indicator;
  real <lower = 0> event_time;

  // phi_{2 \cap 3} -- from the longitudinal submodel
  // used to calculate m'_{i}(t)
  real <lower = 0> breakpoint;
  vector [n_segments] eta_slope; // first segment is eta^{b} -- eta 'before'
}

model {
  if (event_indicator == 1) { // add the hazard
    target += log(hazard_gamma) + (hazard_gamma - 1) * log(event_time);

    if (event_time < breakpoint) {
      target += alpha * eta_slope[1];
    } else {
      target += alpha * eta_slope[2];
    }
  }

  if (event_time < breakpoint) {
    real temp_lower = exp(alpha * eta_slope[1]) * pow(event_time, hazard_gamma);
    target += temp_lower;
  } else {
    real temp_upper = exp(alpha * eta_slope[1]) * pow(breakpoint, hazard_gamma) +
      exp(alpha * eta_slope[2]) * (pow(event_time, hazard_gamma) - pow(breakpoint, hazard_gamma));
    target += temp_upper;
  }
}
