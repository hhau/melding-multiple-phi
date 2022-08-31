data {
 vector [2] lower_boundary;
 vector [2] upper_boundary;
}

parameters {
  vector <lower = lower_boundary, upper = upper_boundary> [2] x;
}

model {
  target += normal_lpdf(x | 0, 1);
}