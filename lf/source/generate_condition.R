#%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Generate Condition ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(latentFactoR); library(mvtnorm); library(compiler)

# Recoverable
recoverable <- function(simulated)
{

  # Update simulated (if necessary)
  if(is(simulated, "lf_pe")){
    simulated <- list(
      population_correlation = simulated$population_correlation,
      parameters = simulated$original_results$parameters,
      overlap = simulated$overlap
    )
  }

  # Estimate EFA on population correlation matrix
  EFA <- EGAnet:::silent_call(efa(
    data = simulated$population_correlation,
    nfactors = simulated$parameters$factors,
    corr = "pearson"
  ))

  # Set up true loadings
  true_loadings <- simulated$parameters$loadings

  # Align with true loadings
  aligned <- EGAnet:::faAlign_fungible(
    F1 = true_loadings,
    F2 = EFA$loadings,
    Phi2 = EFA$Phi
  )

  # Put back into EFA
  EFA$loadings <- aligned$F2
  EFA$Phi <- aligned$Phi2

  # Set true membership
  true_membership <- abs(true_loadings) >= 0.40

  # Return recoverable
  return(
    omega_index(
      obtain_memberships(aligned$F2, true_membership, simulated$overlap),
      true_membership
    ) == 1
  )

}

# Compile function
recoverable <- cmpfun(recoverable)

# Convergence Criterion
convergence <- function(simulated)
{

  # Set up logicals
  max_residual <- simulated$population_error$max_error < 0.10

  # Unidimensional: CFI is limited by the few degrees of freedom
  # so only the residual and discrepancy criteria are used
  if(simulated$original_results$parameters$factors == 1){
    return(max_residual)
  }

  # Multidimensional
  cfi <- abs(simulated$population_error$fit[["CFI"]] - 0.97) < 0.01

  # Return convergence
  return(all(cfi, max_residual))

}

# Compile function
convergence <- cmpfun(convergence)

# Generate condition
generate_condition <- function(condition)
{

  # Total variables
  total_variables <- condition$NVAR * condition$NFAC

  # Set factor sequence
  factor_sequence <- seq_len(condition$NFAC)

  # Number of overlapping cross-loadings to add
  cross_add <- switch(as.character(condition$NVAR), "4" = 1, "8" = 2, "12" = 2)

  # Set unidimensional flag
  unidimensional_flag <- condition$NFAC == 1

  # Initialize regeneration counts
  count_no_overlap <- 0; count_no_overlap_error <- 0
  count_overlap <- 0; count_overlap_error <- 0

  # Loop over communalities
  while(TRUE){

    # Ensure base is recoverable
    while(TRUE){

      # Generate data
      no_overlap <- try(simulate_factors(
        factors = condition$NFAC, variables = condition$NVAR,
        loadings = condition$LOAD, cross_loadings = 0.00,
        correlations = condition$CORF, sample_size = 1
      ), silent = TRUE)

      # Set error flag
      error_flag <- is(no_overlap, "try-error")

      # Add overlap
      no_overlap$overlap <- FALSE

      # Check for unidimensional
      if(unidimensional_flag){
        break
      }else{

        # Check if base is recoverable
        if(error_flag || !recoverable(no_overlap)){

          # Increase count
          count_no_overlap <- count_no_overlap + 1
          next

        }

        # Set loadings
        loadings <- no_overlap$parameters$loadings

        # Add overlapping
        for(k in factor_sequence){

          # Set row index
          row_index <- seq_len(condition$NVAR) + condition$NVAR * (k - 1)
          row_index <- row_index[order(loadings[row_index, k])[seq_len(cross_add)]]

          # k's own higher-order group, and k's position within that group
          h_source   <- ceiling(k / condition$NFAC)
          pos_source <- ((k - 1) %% condition$NFAC) + 1

          # Cyclic successor position within k's own group only
          pos_target <- (pos_source %% condition$NFAC) + 1
          target <- ((h_source - 1) * condition$NFAC) + pos_target

          # Every item flagged for cross-loading on factor k gets the same
          # single deterministic target
          column_index <- rep(target, cross_add)

          # Loop over loadings to add
          for(l in seq_len(cross_add)){

            # Assign cross-loading magnitude
            loadings[row_index[l], c(k, column_index[l])] <- sort(
              EGAnet:::runif_xoshiro(2, min = 0.40, max = 0.50),
              decreasing = TRUE
            )

          }

        }

        # Generate data
        overlap <- try(simulate_factors(
          factors = condition$NFAC, variables = condition$NVAR,
          loadings = loadings, cross_loadings = 0.00,
          correlations = condition$CORF, sample_size = 1
        ), silent = TRUE)

        # Set error flag
        error_flag <- is(overlap, "try-error")

        # Add overlap
        overlap$overlap <- TRUE

        # Check if base is recoverable
        if(!error_flag && recoverable(overlap)){
          break
        }

        # Increase count
        count_overlap <- count_overlap + 1

      }

    }

    # Set CFI target
    target <- 0.97

    # Set CFI target
    if(unidimensional_flag){

      # Loop by increasing CFI target
      while(TRUE){

        # Target check
        target_flag <- abs(target - 1) < 1e-06

        # Obtain error data
        no_overlap_error <- try(add_population_error(
          lf_object = no_overlap, method = "tkl",
          metric = "cfi", fit_value = target,
          control = list(decay = 0.10)
        ), silent = TRUE)

        # Set error flag
        error_flag <- is(no_overlap_error, "try-error")

        # Add overlap
        no_overlap_error$overlap <- FALSE

        # Check if base is recoverable
        if(
          target_flag ||
          !error_flag &&
          convergence(no_overlap_error)
        ){
          break
        }

        # Increase count
        count_no_overlap_error <- count_no_overlap_error + 1

        # Increase target
        target <- sqrt(target)

      }

      # Check for outer loop break
      if(target_flag) next else break

    }else{

      # Obtain error data
      no_overlap_error <- try(add_population_error(
        lf_object = no_overlap, method = "tkl",
        metric = "cfi", fit_value = target,
        control = list(decay = 0.10)
      ), silent = TRUE)

      # Set error flag
      error_flag <- is(no_overlap_error, "try-error")

      # Add overlap
      no_overlap_error$overlap <- FALSE

      # Check if base is recoverable
      if(
        error_flag ||
        !recoverable(no_overlap_error) ||
        !convergence(no_overlap_error)
      ){

        # Increase count
        count_no_overlap_error <- count_no_overlap_error + 1
        next

      }

      # Add error
      overlap_error <- try(add_population_error(
        lf_object = overlap, method = "tkl",
        metric = "cfi", fit_value = target,
        leave_cross_loadings = TRUE,
        control = list(decay = 0.10)
      ), silent = TRUE)

      # Set error flag
      error_flag <- is(overlap_error, "try-error")

      # Add overlap
      overlap_error$overlap <- TRUE

      # Check if base is recoverable
      if(
        !error_flag &&
        recoverable(overlap_error) &&
        convergence(overlap_error)
      ){
        break
      }

      # Increase count
      count_overlap_error <- count_overlap_error + 1

    }

  }

  # Set skew
  skewness <- latentFactoR:::runif_xoshiro(total_variables, min = 0, max = 2)

  # Obtain Cholesky decomposition
  cholesky_no_overlap <- chol(no_overlap$population_correlation)
  cholesky_no_overlap_error <- chol(no_overlap_error$population_correlation)

  # Check for condition
  if(!unidimensional_flag){
    cholesky_overlap <- chol(overlap$population_correlation)
    cholesky_overlap_error <- chol(overlap_error$population_correlation)
  }

  # Obtain means and Sigma
  means <- rep(0, total_variables); Sigma <- diag(total_variables)

  # Set sample sizes
  N <- c(500, 1000, 2500, 5000)

  # Generate samples
  samples <- lapply(N, function(n){

    # Generate data
    data <- rmvnorm(n = n, mean = means, sigma = Sigma)

    # Create correlated data
    data_no_overlap <- data %*% cholesky_no_overlap
    data_no_overlap_error <- data %*% cholesky_no_overlap_error

    # Check for condition
    if(!unidimensional_flag){
      data_overlap <- data %*% cholesky_overlap
      data_overlap_error <- data %*% cholesky_overlap_error
    }

    # Add skew to the data
    for(v in seq_len(total_variables)){

      # Add skew
      data_no_overlap[,v] <- categorize(
        data_no_overlap[,v], categories = 5, skew_value = skewness[v]
      )
      data_no_overlap_error[,v] <- categorize(
        data_no_overlap_error[,v], categories = 5, skew_value = skewness[v]
      )

      # Check for condition
      if(!unidimensional_flag){
        data_overlap[,v] <- categorize(
          data_overlap[,v], categories = 5, skew_value = skewness[v]
        )
        data_overlap_error[,v] <- categorize(
          data_overlap_error[,v], categories = 5, skew_value = skewness[v]
        )
      }

    }

    # Set up samples
    result <- list(no_overlap = data_no_overlap, no_overlap_error = data_no_overlap_error)

    # Check for unidimensional
    if(!unidimensional_flag){
      result$overlap <- data_overlap
      result$overlap_error <- data_overlap_error
    }

    # Return samples
    return(result)

  })

  # Set up final list
  result_list <- list()

  # Add samples to conditions
  no_overlap$data <- structure(lapply(samples, function(x){x$no_overlap}), names = N)
  no_overlap$regenerate <- count_no_overlap
  result_list$no_overlap <- no_overlap
  rm(no_overlap)
  no_overlap_error$data <- structure(lapply(samples, function(x){x$no_overlap_error}), names = N)
  no_overlap_error$regenerate <- count_no_overlap_error
  result_list$no_overlap_error <- no_overlap_error
  rm(no_overlap_error)

  # Check for unidimensional
  if(!unidimensional_flag){
    overlap$data <- structure(lapply(samples, function(x){x$overlap}), names = N)
    overlap$regenerate <- count_overlap
    result_list$overlap <- overlap
    rm(overlap)
    overlap_error$data <- structure(lapply(samples, function(x){x$overlap_error}), names = N)
    overlap_error$regenerate <- count_overlap_error
    result_list$overlap_error <- overlap_error
    rm(overlap_error)
  }



  # Return final list
  return(result_list)

}

# Compile function
generate_condition <- cmpfun(generate_condition)
