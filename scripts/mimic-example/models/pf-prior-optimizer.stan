data {
  int <lower = 1> n_prior_samples;
  int <lower = 1, upper = 2> event_indicator [n_prior_samples];
  vector <lower = 0> [n_prior_samples] event_time;
  real <lower = 0> length_of_stay;
}

transformed data {
  vector <lower = 0, upper = 1> [n_prior_samples] event_time_scaled;
  event_time_scaled = event_time / length_of_stay;
}

parameters {
  real <lower = 0, upper = 1> weight;
  real <lower = 0> beta_alpha;
  real <lower = 0> beta_beta;
}

model {
  for (ii in 1 : n_prior_samples) {
    if (event_indicator[ii] == 1) {
      target += beta_lpdf(event_time_scaled[ii] | beta_alpha, beta_beta);
      target += log(weight);
    } else if (event_indicator[ii] == 2) {
      target += log(1 - weight);
    }
    // this jacobian is constant and so technically could be ignored, but i'd
    // forget that i was being clever and would mess something up down the road
    target += -log(length_of_stay);
  }
}
