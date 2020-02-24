data {
  matrix [2, 2] sigma_mat;
}

parameters {
  real phi_12;
  real phi_23;
}

transformed parameters {
  vector [2] phi = [phi_12, phi_23]';
}

model {
  target += normal_lpdf(phi_12 | -2.0, 1.0);
  target += multi_normal_lpdf(phi | [0.0, 0.0], sigma_mat);
  target += normal_lpdf(phi_23 | 2.0, 1.0);
}
