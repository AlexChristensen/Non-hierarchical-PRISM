#%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Generate Condition ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(psychNetwoRk); library(sbm); library(igraph)
library(mvtnorm); library(psych); library(GPArotation); library(compiler)

# Recoverable
recoverable <- function(simulated)
{

  # Check for pseudo Omega
  if(is(simulated, "pn_pe")){
    Omega <- simulated$population_error$Omega_pseudo
    blocks <- simulated$original_results$parameters$blocks
    membership <- simulated$original_results$parameters$membership
  }else{
    Omega <- simulated$parameters$network
    blocks <- simulated$parameters$blocks
    membership <- simulated$parameters$membership
  }

  # Ensure symmetric
  Omega <- (Omega + t(Omega)) / 2

  # Estimate SBM
  SBM <- estimateSimpleSBM(
    netMat = Omega, model = "gaussian",
    directed = FALSE, estimOptions = list(
      nbCores = 1, verbosity = 0, plot = FALSE,
      nbBlocksRange = blocks
    )
  )

  # Ensure population communities
  block_index <- which(SBM$storedModels$nbBlocks == blocks)

  # Check the blocks exists
  if(length(block_index) == 0){
    return(FALSE)
  }

  # Update model
  SBM$setModel(block_index)

  # Return recoverable
  return(abs(compare(SBM$memberships, membership, method = "adjusted.rand") - 1) < 1e-06)

}

# Compile function
recoverable <- cmpfun(recoverable)

# Convergence Criterion
convergence <- function(simulated)
{

  # Set up logicals
  max_residual <- simulated$population_error$max_error < 0.10
  max_discrepancy <- max(abs(simulated$population_error$discrepancy)) < 0.10

  # Unidimensional: CFI is limited by the few degrees of freedom
  # so only the residual and discrepancy criteria are used
  if(simulated$original_results$parameters$blocks == 1){
    return(all(max_residual, max_discrepancy))
  }

  # Multidimensional
  cfi <- abs(simulated$population_error$fit[["CFI"]] - 0.97) < 0.01

  # Return convergence
  return(all(cfi, max_residual, max_discrepancy))

}

# Compile function
convergence <- cmpfun(convergence)

# EFA
efa_fit <- function(R, nfactors, randomStarts = 100L, ...)
{

  # Estimate EFA on population correlation matrix
  EFA <- EGAnet:::silent_call(fa(
    r = R, n.obs = 8e09, nfactors = nfactors,
    rotate = "none", fm = "ml"
  ))

  # Check for whether rotation is necessary
  if(nfactors != 1){

    # Rotate loadings
    rotated <- GPArotation::geominQ(
      EFA$loadings, randomStarts = randomStarts, eps = switch(
        as.character(nfactors),
        "2" = 0.0001,
        "3" = 0.001,
        0.01
      )
    )

    # Set up correlations
    implied_R <- rotated$loadings %*% tcrossprod(rotated$Phi, rotated$loadings)

  }else{

    # Set loadings
    loadings <- EFA$loadings[,seq_len(nfactors), drop = FALSE]

    # Set up correlations
    implied_R <- tcrossprod(loadings)

  }

  # Set diagonal to one
  diag(implied_R) <- 1

  # Return
  return(implied_R)

}

# Compile function
efa_fit <- cmpfun(efa_fit)

# Generate condition
generate_condition <- function(condition)
{

  # Set density matrix
  DM <- matrix(condition$BETWEEN, nrow = condition$BLOCKS, ncol = condition$BLOCKS)
  diag(DM) <- condition$WITHIN

  # Set unidimensional flag
  unidimensional_flag <- condition$BLOCKS == 1

  # Initialize regeneration counts
  count_base <- 0; count_error <- 0

  # Generate until recoverable
  while(TRUE){

    # Generate until recoverable
    while(TRUE){

      # Generate data
      simulated_base <- try(simulate_sbm(
        nodes = condition$NODES, blocks = condition$BLOCKS,
        density_matrix = DM, diffusion = condition$DIFFUSION,
        sample_size = 1, negative_proportion = 0
      ), silent = TRUE)

      # Check if base is recoverable
      if(!is(simulated_base, "try-error")){

        # Check for unidimensional or recoverable
        if(unidimensional_flag || recoverable(simulated_base)){
          break
        }

      }

      # Increase count
      count_base <- count_base + 1

    }

    # Set CFI target
    target <- 0.97

    # Set CFI target
    if(unidimensional_flag){

      # Loop by increasing CFI target
      while(TRUE){

        # Check for target of one
        if(abs(target - 1) < 1e-06){
          next
        }

        # Obtain error data
        simulated_error <- try(add_population_error(
          pn_object = simulated_base, method = "tkl", metric = "cfi",
          fit_value = target, control = list(decay = 0.10)
        ), silent = TRUE)

        # Check if base is recoverable
        if(!is(simulated_error, "try-error") && convergence(simulated_error)){
          break
        }

        # Increase count
        count_error <- count_error + 1

        # Increase target
        target <- sqrt(target)

      }

    }else{

      # Obtain error data
      simulated_error <- try(add_population_error(
        pn_object = simulated_base, method = "tkl", metric = "cfi",
        fit_value = target, control = list(decay = 0.10)
      ), silent = TRUE)

      # Check if base is recoverable
      if(
        !is(simulated_error, "try-error") &&
        convergence(simulated_error) &&
        recoverable(simulated_error)
      ){
        break
      }

      # Increase count
      count_error <- count_error + 1

    }

  }

  # Obtain EFA fits
  ## Base
  efa_base <- efa_fit(
    R = simulated_base$population_correlation,
    nfactors = simulated_base$parameters$blocks
  )
  ## Error
  efa_error <- efa_fit(
    R = simulated_error$population_correlation,
    nfactors = simulated_base$parameters$blocks
  )

  # Set degrees of freedom
  total_variables <- sum(simulated_base$parameters$nodes)
  df <- (total_variables * (total_variables - 1) / 2) - (total_variables * condition$BLOCKS)

  # Obtain fit indices
  ## Base
  efa_base <- latentFactoR::obtain_fit_indices(
    R_model = efa_base,
    R = simulated_base$population_correlation,
    sample_size = Inf,
    df = df
  )
  ## Base
  efa_error <- latentFactoR::obtain_fit_indices(
    R_model = efa_error,
    R = simulated_base$population_correlation,
    sample_size = Inf,
    df = df
  )

  # Obtain number of variables
  n_variables <- length(simulated_base$parameters$membership)

  # Set skew
  skewness <- psychNetwoRk:::runif_xoshiro(n_variables, min = 0, max = 2)

  # Obtain Cholesky decomposition
  cholesky_base <- chol(simulated_base$population_correlation)
  cholesky_error <- chol(simulated_error$population_correlation)

  # Obtain means and Sigma
  means <- rep(0, n_variables); Sigma <- diag(n_variables)

  # Set sample sizes
  N <- c(500, 1000, 2500, 5000)

  # Generate samples
  samples <- lapply(N, function(n){

    # Generate data
    data <- rmvnorm(n = n, mean = means, sigma = Sigma)

    # Create correlated data
    data_base <- data %*% cholesky_base
    data_error <- data %*% cholesky_error

    # Add skew to the data
    for(v in seq_len(n_variables)){

      # Add skew
      data_base[,v] <- categorize(data_base[,v], categories = 5, skew_value = skewness[v])
      data_error[,v] <- categorize(data_error[,v], categories = 5, skew_value = skewness[v])

    }

    # Return samples
    return(list(base = data_base, error = data_error))

  })

  # Set up final list
  result_list <- list()

  # Add samples to conditions
  simulated_base$data <- structure(lapply(samples, function(x){x$base}), names = N)
  simulated_base$regenerate <- count_base
  simulated_base$efa_fit <- efa_base
  result_list$base <- simulated_base
  rm(simulated_base)
  simulated_error$data <- structure(lapply(samples, function(x){x$error}), names = N)
  simulated_error$regenerate <- count_error
  simulated_error$efa_fit <- efa_error
  result_list$error <- simulated_error
  rm(simulated_error)

  # Return final list
  return(result_list)

}

# Compile function
generate_condition <- cmpfun(generate_condition)
