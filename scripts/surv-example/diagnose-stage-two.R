library(bayesplot)
library(ggplot2)

source("scripts/common/plot-settings.R")

phi_samples <- readRDS("rds/surv-example/stage-two-phi-samples.rds")

mcmc_trace(phi_samples)
rstan::monitor(phi_samples, warmup = 1) 
# so we have 14 Tail_ESS per 1000 samples, to get to ~500, we need
# (500 / 14) * 1000 ~~ 35,000, and we can chuck on 
mcmc_acf(phi_samples)
