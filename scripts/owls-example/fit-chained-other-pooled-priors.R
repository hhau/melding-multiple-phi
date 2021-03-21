# Because of the flat / uninformative nature of the priors, we can sample from
# the chained melded posterior that uses either pooling by resampling
# the PoE posterior under the new pooled prior. It will have almost no effect
# but will at least be a complete set of melded posteriors.
library(abind)
library(parallel)

source("scripts/common/mcmc-util.R")
source("scripts/common/logger-setup.R")

flog.info(
  "owls-example-other-pooled-priors: loading data and setup",
  name = base_filename 
)

melded_posterior_samples <- readRDS(
  "rds/owls-example/melded-posterior-samples.rds"
)

n_chain <- 6
n_iter_per_chain <- 2e4

n_prev_chain <- dim(melded_posterior_samples)[2]
n_prev_iter <- dim(melded_posterior_samples)[1]

phi_12_names <- c("v[1]", "v[2]")
phi_23_names <- c("fec")
phi_names <- c(phi_12_names, phi_23_names)
n_param <- length(phi_names)

# pooling weights and targets.
lambda_log_pooling <- rep(1 / 2, times = 3)
lambda_linear_pooling <- rep(1 / 2, times = 4)

q_pool_log <- function(phi_vec, log = TRUE) {
  log_res <- sum(dnorm(phi_vec[1 : 2], sd = 100, log = TRUE)) *
    sum(lambda_log_pooling[1 : 2])

  log_res <- log_res + dunif(phi_vec[3], min = 0, max = 10, log = TRUE) *
    sum(lambda_log_pooling[2 : 3])

  if (log) return(log_res) else return(exp(log_res))
}

target_log_pooling <- function(phi_vec, log = TRUE) {
  log_res <- q_pool_log(phi_vec)
  log_res <- log_res - 2 * sum(dnorm(phi_vec[1 : 2], sd = 100, log = TRUE))
  log_res <- log_res - 2 * dunif(phi_vec[3], min = 0, max = 10, log = TRUE)

  log_res <- unname(log_res)
  if (log) return(log_res) else return(exp(log_res))
}

# the uniform prior is irrelevant in the linear pooling case, actually I guess
# it was also irrelevant in the log pooling case as well.
# in this example, linear pooling is completely independent of the pooling
# weights, which is kinda interesting.
# with log pooling weights lambda_1 = lambda_2 = 1 / 2, this is exactly the 
# same as log pooling
q_pool_lin <- function(phi_vec, log = TRUE) {

}

target_lin_pooling <- function(phi_vec, log = TRUE) {
  log_res <- -sum(dnorm(phi_vec[1 : 2], mean = 0, sd = 100, log = TRUE))
  if (log) log_res else exp(log_res)
}

# start the mclapply
res_log <- mclapply(1 : n_chain, mc.cores = 6, function(chain_id) {
  flog.info(
    sprintf("owls-example-log-pooling: starting chain-%d", chain_id),
    name = base_filename
  )
  
  previous_stage_indices <- array(
    NA,
    dim = c(n_iter_per_chain, 1, 2),
    dimnames = list(
      sprintf("iteration_%d", 1 : n_iter_per_chain),
      sprintf("chain_%d", chain_id),
      c("prev_iteration_index", "prev_chain_index")
    )
  )
  
  phi_samples <- array(
    NA,
    dim = c(n_iter_per_chain, 1, n_param),
    dimnames = list(
      sprintf("iteration_%d", 1 : n_iter_per_chain),
      sprintf("chain_%d", chain_id),
      phi_names
    )
  )

  iter_init <- sample.int(n = n_prev_iter, size = 1)
  chain_init <- sample.int(n = n_prev_chain, size = 1)
  previous_stage_indices[1, 1, ] <- c(iter_init, chain_init)
  
  phi_samples[1, 1, phi_names] <- melded_posterior_samples[
    iter_init,
    chain_init,
    phi_names
  ]
  
  for (ii in 2 : n_iter_per_chain) {
    # start the inner loop
    iter_prop <- sample.int(n = n_prev_iter, size = 1)
    chain_prop <- sample.int(n = n_prev_chain, size = 1)
    iter_curr <- previous_stage_indices[ii - 1, 1, "prev_iteration_index"]
    chain_curr <- previous_stage_indices[ii - 1, 1, "prev_chain_index"]
    
    phi_prop <- melded_posterior_samples[iter_prop, chain_prop, phi_names]
    phi_curr <- phi_samples[ii - 1, 1, phi_names]
    log_alpha <- target_log_pooling(phi_prop) - target_log_pooling(phi_curr)
    
    if (runif(1) < exp(log_alpha)) {
      previous_stage_indices[ii, 1, ] <- c(iter_prop, chain_prop)
      phi_samples[ii, 1, ] <- phi_prop
    } else {
      previous_stage_indices[ii, 1, ] <- previous_stage_indices[ii - 1, 1, ]
      phi_samples[ii, 1, ] <- phi_curr
    }
    
    if (ii %% 500 == 0) {
      flog.info(
        sprintf(
          "owls-example-log-pooling: chain-%d -- iteration-%d", 
          chain_id, 
          ii
        ),
        name = base_filename
      )
    }
  }

  res_list <- list(
    previous_stage_indices = previous_stage_indices,
    phi_samples = phi_samples
  )

  return(res_list)
})

res_lin <- mclapply(1 : n_chain, mc.cores = 6, function(chain_id) {
  flog.info(
    sprintf("owls-example-lin-pooling: starting chain-%d", chain_id),
    name = base_filename
  )
  
  previous_stage_indices <- array(
    NA,
    dim = c(n_iter_per_chain, 1, 2),
    dimnames = list(
      sprintf("iteration_%d", 1 : n_iter_per_chain),
      sprintf("chain_%d", chain_id),
      c("prev_iteration_index", "prev_chain_index")
    )
  )
  
  phi_samples <- array(
    NA,
    dim = c(n_iter_per_chain, 1, n_param),
    dimnames = list(
      sprintf("iteration_%d", 1 : n_iter_per_chain),
      sprintf("chain_%d", chain_id),
      phi_names
    )
  )
  
  iter_init <- sample.int(n = n_prev_iter, size = 1)
  chain_init <- sample.int(n = n_prev_chain, size = 1)
  previous_stage_indices[1, 1, ] <- c(iter_init, chain_init)
  
  phi_samples[1, 1, phi_names] <- melded_posterior_samples[
    iter_init,
    chain_init,
    phi_names
  ]
  
  for (ii in 2 : n_iter_per_chain) {
    # start the inner loop
    iter_prop <- sample.int(n = n_prev_iter, size = 1)
    chain_prop <- sample.int(n = n_prev_chain, size = 1)
    iter_curr <- previous_stage_indices[ii - 1, 1, "prev_iteration_index"]
    chain_curr <- previous_stage_indices[ii - 1, 1, "prev_chain_index"]
    
    phi_prop <- melded_posterior_samples[iter_prop, chain_prop, phi_names]
    phi_curr <- phi_samples[ii - 1, 1, phi_names]
    log_alpha <- target_lin_pooling(phi_prop) - target_lin_pooling(phi_curr)
    
    if (runif(1) < exp(log_alpha)) {
      previous_stage_indices[ii, 1, ] <- c(iter_prop, chain_prop)
      phi_samples[ii, 1, ] <- phi_prop
    } else {
      previous_stage_indices[ii, 1, ] <- previous_stage_indices[ii - 1, 1, ]
      phi_samples[ii, 1, ] <- phi_curr
    }
    
    if (ii %% 5000 == 0) {
      flog.info(
        sprintf(
          "owls-example-lin-pooling: chain-%d -- iteration-%d", 
          chain_id, 
          ii
        ),
        name = base_filename
      )
    }
  }
  
  res_list <- list(
    previous_stage_indices = previous_stage_indices,
    phi_samples = phi_samples
  )

  return(res_list)
})

bind_named_sublists <- function(outer_list, name) {
  lapply(outer_list, function(sub_list) sub_list[[name]]) %>% 
    abind(along = 2)
}

phi_samples_log_pooling <- bind_named_sublists(res_log, "phi_samples")
indices_log_pooling <- bind_named_sublists(res_log, "previous_stage_indices")

phi_samples_lin_pooling <- bind_named_sublists(res_lin, "phi_samples")
indices_lin_pooling <- bind_named_sublists(res_lin, "previous_stage_indices")

flog.info(
  "owls-example-other-pooled-priors: writing samples to disk",
  name = base_filename
)

saveRDS(
  object = phi_samples_log_pooling,
  file = "rds/owls-example/melded-phi-samples-log-pooling.rds"
)

saveRDS(
  object = indices_log_pooling,
  file = "rds/owls-example/melded-indices-log-pooling.rds"
)

saveRDS(
  object = phi_samples_lin_pooling,
  file = "rds/owls-example/melded-phi-samples-lin-pooling.rds"
)

saveRDS(
  object = indices_lin_pooling,
  file = "rds/owls-example/melded-indices-lin-pooling.rds"
)

# renamer <- function(x) {
#   res <- c("v[1]" = "alpha[0]", "v[2]" = "alpha[2]", "fec" = "rho")
#   res[x]
# }
# 
# plot_worst_pars(phi_samples_log_pooling, facet_name_value_pairs = renamer)
