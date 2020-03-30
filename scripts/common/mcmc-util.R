mcmc_list_to_array <- function(mcmc_list) {
  stopifnot(class(mcmc_list) == "mcmc.list")
  
  n_iter <- nrow(mcmc_list[[1]])
  n_par <- ncol(mcmc_list[[1]])
  n_chain <- length(mcmc_list)
  par_names <- dimnames(mcmc_list[[1]])[[2]]
  
  results_array <- array(
    data = NA,
    dim = c(n_iter, n_chain, n_par),
    dimnames = list(
      NULL,
      sprintf("chain_%d", 1:n_chain),
      par_names
    )
  )
  
  for (chain_index in 1:n_chain) {
    results_array[, chain_index,] <- mcmc_list[[chain_index]]
  }
  
  return(results_array)
  
}

array_to_mcmc_list <- function(array) {
  stopifnot(length(dim(array)) == 3)

  n_chain <- dim(array)[2]
  res <- lapply(1 : n_chain, function(chain_id) {
    coda::mcmc(array[ , chain_id, ])
  })

  return(coda::as.mcmc.list(res))
}