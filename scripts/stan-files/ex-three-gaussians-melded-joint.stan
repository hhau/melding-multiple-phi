data {
  int n_1;
  vector [n_1] y_1;

  int n_2;
  vector [n_2] y_2;

  int n_3;
  vector [n_3] y_3;
}

parameters {
  real phi_12;
  real phi_23;

  real <lower = 0> sigma_y_1;
  real <lower = 0> sigma_y_2;
  real <lower = 0> sigma_y_3;

}

model {
  // model 1
  target += normal_lpdf(y_1 | phi_12, sigma_y_1);
  target += normal_lpdf(phi_12 | 1.0, 1.0);
  target += normal_lpdf(sigma_y_1 | 0.0, 1.0);

  // model 2
  target += normal_lpdf(y_2 | phi_12 + phi_23, sigma_y_2);
  target += normal_lpdf(phi_12 | 1.0, 1.0);
  target += normal_lpdf(phi_23 | 2.0, 1.0);
  target += normal_lpdf(sigma_y_2 | 0.0, 1.0);

  // model 3
  target += normal_lpdf(y_3 | phi_23, sigma_y_3);
  target += normal_lpdf(phi_23 | 1.0, 1.0);
  target += normal_lpdf(sigma_y_3 | 0.0, 1.0);

  // marginalise out all the priors on phis
  target += -1.0 * normal_lpdf(phi_12 | 1.0, 1.0);
  target += -1.0 * normal_lpdf(phi_12 | 1.0, 1.0);
  target += -1.0 * normal_lpdf(phi_23 | 2.0, 1.0);
  target += -1.0 * normal_lpdf(phi_23 | 1.0, 1.0);

  // Pooled prior - equal weighted log pooling 
  target += 0.5 * normal_lpdf(phi_12 | 1.0, 1.0);
  target += 0.5 * normal_lpdf(phi_12 | 1.0, 1.0);

  target += 0.5 * normal_lpdf(phi_23 | 2.0, 1.0);
  target += 0.5 * normal_lpdf(phi_23 | 1.0, 1.0);

}
