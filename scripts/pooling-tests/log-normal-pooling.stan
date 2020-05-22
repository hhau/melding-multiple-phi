data {
  matrix [2, 2] sigma_mat;
  vector [3] lambda_vec;
}

parameters {
  real phi_12;
  real phi_23;
}

transformed parameters {
  vector [2] phi = [phi_12, phi_23]';
}

model {
  target += lambda_vec[1] * normal_lpdf(phi_12 | -2.0, 1.0);
  target += lambda_vec[2] * multi_normal_lpdf(phi | [0.0, 0.0], sigma_mat);
  target += lambda_vec[3] * normal_lpdf(phi_23 | 2.0, 1.0);
}
