data {
  // Number of time points
  int <lower = 0> T;

  // Observations
  int <lower = 0> N_breeding_females [T];
  int <lower = 0> N_offspring [T];
}

parameters {
  // In other variants of this model, this is also a function of t.
  // Here, we take the model with the highest evidence 
  // from Finke which has it as constant.
  real <lower = 0, upper = 10> rho;
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
}