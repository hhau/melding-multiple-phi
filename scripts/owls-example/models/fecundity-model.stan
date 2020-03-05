data {
  // Number of time points
  int <lower = 0> T;

  // Observations
  int <lower = 0> N_breeding_females [T];
  int <lower = 0> N_offspring [T];
}

parameters {
  // In other variants of this model, this is also a function of t.
  // Here, we take the highest evidence model from Finke and make it constant.
  real <lower = 0> rho;
}

transformed parameters {
  real <lower = 0> product [T];

  for (tt in 1 : T) {
    product[tt] = N_breeding_females[tt] * rho;
  }
}

model {
  // Model as per Finke (2019)
  target += poisson_lpmf(N_offspring | product);
  // Original paper uses uniform[0, 10] prior, but that is silly.
  // Here's a more sensible half normal that still puts most of its mass between
  // 0 and 10 
  target += normal_lpdf(rho | 0.0, 3.0);
}