functions {
  vector pow_vec(vector x, real y) {
    int N = rows(x);
    vector [N] res;
    for (nn in 1 : N) {
      res[nn] = pow(x[nn], y);
    }
    return res;
  }

  real integrand_stable (real x, real xc, real[] theta, real[] x_r, int[] x_i) {
    real hazard_gamma = theta[1];
    real alpha = theta[2];
    real beta_one = theta[3];
    real sign = (alpha / fabs(alpha)) * (beta_one / fabs(beta_one));
    real lf = log(fabs(alpha)) + log(fabs(beta_one)) - log(hazard_gamma) +  
      hazard_gamma * log(x) + alpha * beta_one * x;
    real res = sign * exp(lf);
    
    if(is_nan(res)) {
      print(" x: ", x, " beta_1:", beta_one, " f(x):", res, " logf(x):", lf);
    }

    return(res);
  }
}

data {
  // global constant
  int <lower = 1> n_patients;

  // Y_{2} - baseline data
  vector [n_patients] baseline_measurement;
  real log_crude_event_rate;

  // phi_{1 \cap 2} -- from the event time submodel
  // assuming no censoring for the moment
  vector <lower = 0> [n_patients] event_time;
  
  // phi_{2 \cap 3} -- from the longitudinal submodel
  // used to calculate m_{i}(t)
  matrix [n_patients, 2] long_beta;
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  // intercept
  real beta_zero;

  // baseline coefficient
  real beta_one;

  // baseline hazard parameter(s) (Weibull gamma)
  real <lower = 0> hazard_gamma;

  // longitudinal associative strength (alpha)
  real alpha;
}

transformed parameters {
  vector [n_patients] z_common;
  vector [n_patients] surv_prob_term_one;
  vector [n_patients] surv_prob_term_two;

  z_common = beta_zero + 
    baseline_measurement * beta_one + 
    alpha * long_beta[, 1];

  surv_prob_term_one = hazard_gamma * exp(z_common);

  for (ii in 1 : n_patients) {
    real I2; // temp stuff to store intermediate
    real w;
    I2 = integrate_1d(
      integrand_stable,
      0.0,
      event_time[ii],
      {hazard_gamma, alpha, long_beta[ii, 2]},
      x_r,
      x_i,
      0.01
    );

    w = exp(
      hazard_gamma * log(event_time[ii]) + 
      alpha * long_beta[ii, 2] * event_time[ii] - 
      log(hazard_gamma)
    );

    surv_prob_term_two[ii] = w - I2;
  }
}

model {
  // Weibull model, no censoring.
  // log hazard
  target += n_patients * log(hazard_gamma);
  target += (hazard_gamma - 1) * sum(log(event_time));
  target += sum(z_common) + alpha * sum(long_beta[, 2] .* event_time);

  // log survival probability
  target += -sum(surv_prob_term_one .* surv_prob_term_two);

  // priors 
  target += normal_lpdf(beta_zero | log_crude_event_rate, 1);
  target += normal_lpdf(beta_one | 0, 1);
  target += normal_lpdf(hazard_gamma | 0, 1);
  target += normal_lpdf(alpha | 0, 2);
}
