// this model is only for updating an individual at a time in a very specific
// MH-within-Gibbs scheme, it does not compute the full log-posterior.
// only \phi_{1 \cap 2, i} \mid \phi_{2 \cap 3}, \psi_{2}, Y_{2}
// very much up to proportionality -- many constants are dropped.

functions {
  real integrand_stable (real x, real xc, real[] theta, real[] x_r, int[] x_i) {
    real hazard_gamma = theta[1];
    real alpha = theta[2];
    real beta_one = theta[3];
    real sign = (alpha / fabs(alpha)) * (beta_one / fabs(beta_one));
    real lf = log(fabs(alpha)) + log(fabs(beta_one)) - log(hazard_gamma) +  
      hazard_gamma * log(x) + alpha * beta_one * x;
    real res = sign * exp(lf);
    
    return(res);
  }
}

data {
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  // baseline hazard parameter(s) (Weibull gamma)
  real <lower = 0> hazard_gamma;

  // longitudinal associative strength (alpha)
  real alpha;

  // phi_{1 \cap 2} -- from the event time submodel
  // note that event_indicator will be a vector of ones n zeros, but 
  // we can't put the bounds on it because they will have log_prob of -inf
  real <lower = 0> event_time;
  real event_indicator;
  real eta_one;
}

transformed parameters {
  real I;
  real I2; // temp stuff to store intermediate
  real w;
  I2 = integrate_1d(
    integrand_stable,
    0.0,
    event_time,
    {hazard_gamma, alpha, eta_one},
    x_r,
    x_i,
    0.0001
  );  
  
  w = exp(
    hazard_gamma * log(event_time) + 
    alpha * eta_one * event_time - 
    log(hazard_gamma)
  );
  // note that various constants have been pulled into I2 to improve 
  // the numerics, which is why this doesn't appear to match the maths
  I = w - I2;
}

model {
  // Weibull model log hazard
  // note that this model is never used to compute a gradient
  // which is why I can do this branching 
  if (event_indicator == 1) {
    target += (hazard_gamma - 1) * log(event_time);
    target += alpha * eta_one * event_time;
  }
  
  // log survival probability
  target += -I;
}
