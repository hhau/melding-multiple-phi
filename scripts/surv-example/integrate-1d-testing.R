suppressMessages(library(cmdstanr))

cmdstan_model <- cmdstan_model(
  "/Users/amanderson/OneDrive - The Alan Turing Institute/PhD/Projects/2019-10-11_multiple-phi/scripts/surv-example/models/submodel-two-phi-step.stan"
)

stan_data <- list(
  n_patients = 16L, 
  baseline_measurement = c(
   -0.935695521120025, -2.12761603089136, -0.876325387003469, 0.931350429946891,
    0.422731897638725, -1.30375210258178,  1.85158267797611, -0.289448169455647, 
    0.325903831278964, -1.2757407765375,   1.58794213975781,  1.59307842403321, 
   -0.0584097720952887, 1.80357260097634, -2.30428952138529,  0.655115279462301
 ), 
 log_crude_event_rate = 3.17644992456411
)

sink(file = "logs/numerics.log", type = "output")
sink(file = "logs/numerics.log", type = "message")
res <- cmdstan_model$sample(
  data = stan_data,
  chains = 1,
  iter_warmup = 1000,
  iter_sampling = 200,
  init = 0.01,
  adapt_delta = 0.99
)
sink(type = "message")
sink()     

# some prefab examples:

# Example - decaying behaviour
event_time <- 0.333964
hazard_gamma <- 0.280632
alpha <- 1.71436
beta_one <- 1.641

# Chain 1 Rejecting initial value:
# Chain 1   Error evaluating the log probability at the initial value.
# Chain 1 Exception: integrate: error estimate of integral 0.0568744 exceeds the
# given relative tolerance times norm of integral

f_integrand <- function(x) {
  x^(hazard_gamma - 1) * exp(alpha * beta_one * x)
}

f_alt <- function(x) {
  exp((hazard_gamma - 1) * log(x) + alpha * beta_one * x)
}

curve(f_integrand, from = 0, to = event_time)
curve(f_alt, from = 0, to = event_time, add = T, col = "blue")

integrate(f_integrand, lower = 0, upper = event_time)
integrate(f_alt, lower = 0, upper = event_time)

# And another - increasing behaviour
event_time <- 0.141916
hazard_gamma <- 4.43969
alpha <- -1.20662 
beta_one <- -0.81199 

#Chain 1 Exception: integrate: error estimate of integral 7.81929e-08 exceeds

f_integrand <- function(x) {
  x^(hazard_gamma - 1) * exp(alpha * beta_one * x)
}

f_alt <- function(x) {
  exp((hazard_gamma - 1) * log(x) + alpha * beta_one * x)
}

curve(f_integrand, from = 0, to = event_time)
curve(f_alt, from = 0, to = event_time, add = T, col = "blue")

integrate(f_integrand, lower = 0, upper = event_time)
integrate(f_alt, lower = 0, upper = event_time)
