data {
  int <lower = 1> n_patients;

  // event times and censoring indicators
  int <lower = 0, upper = 1> event_indicator [n_patients];
  // we are assuming that event_time is censored at t = 1
  real <lower = 0, upper = 1> event_time [n_patients];
  
  // trajectory model parameters
  // we may want to change the type here (vector <-> real array)
  // also at the moment this is only for one time varying covariate
  // trajectory. Matrix type / n_covariates also required if we want to 
  // do more
  vector [n_patients] traj_beta_zero;
  vector [n_patients] traj_beta_one;
}

transformed data {
  int <lower = 0> n_covariates = 1;
}

parameters {
  // coefficients of the trajectories (inside the survival model)
  // there is also an intercept here?
  vector [n_covariates] traj_alpha;

  // baseline hazard parameter(s) (weibull lambda)
  vector <lower = 0> hazard_lambda;

}

model {
  // ether way we have to add the survival prob to the log posterior

  // if the observation is uncensored then we also have to add the harzard

  // we don't have to do the time varing/trajectory model here? assuming it is
  // done in it's own submodel
  // two-stage joint modelling, with uncertain event times. 

}