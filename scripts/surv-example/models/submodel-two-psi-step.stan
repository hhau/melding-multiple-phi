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
    
    return(res);
  }
}

data {
  // global constant
  int <lower = 1> n_patients;
  int <lower = 1> n_long_beta;

  // Y_{2} - baseline data
  vector [n_patients] baseline_measurement;
  real log_crude_event_rate;

  // phi_{1 \cap 2} -- from the event time submodel
  int <lower = 0, upper = 1> event_indicator [n_patients]; 
  vector <lower = 0> [n_patients] event_time;
  
  // phi_{2 \cap 3} -- from the longitudinal submodel
  // used to calculate m_{i}(t)
  matrix [n_patients, n_long_beta] long_beta;
}

transformed data {
  int <lower = 0, upper = n_patients> n_event_patients = sum(event_indicator);
  real x_r[0];
  int x_i[0];
}

parameters {
  // intercept
  real theta_zero;

  // baseline coefficient
  real theta_one;

  // baseline hazard parameter(s) (Weibull gamma)
  real <lower = 0> hazard_gamma;

  // longitudinal associative strength (alpha)
  real alpha;
}

transformed parameters {
  vector [n_patients] z_common;
  vector [n_patients] surv_prob_term_one;
  vector [n_patients] surv_prob_term_two;

  z_common = theta_zero + 
    baseline_measurement * theta_one + 
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
      0.0001
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
  // Weibull model log hazard
  // note that this model is never used to compute a gradient
  // which is why I can do this branching 
  for (ii in 1 : n_patients) {
    if (event_indicator[ii] == 1) {
      target += hazard_gamma;
      target += (hazard_gamma - 1) * log(event_time[ii]);
      target += z_common[ii];
      target += alpha * long_beta[ii, 2] * event_time[ii];
    }
  }

  // log survival probability
  target += -sum(
    surv_prob_term_one .* surv_prob_term_two
  );

  // priors 
  target += normal_lpdf(theta_zero | log_crude_event_rate, 1);
  target += normal_lpdf(theta_one | 0, 1);
  target += normal_lpdf(hazard_gamma | 0, 1);
  target += normal_lpdf(alpha | 0, 1);
}
