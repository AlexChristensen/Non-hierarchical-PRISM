#%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Generate Condition ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(latentFactoR); library(mvtnorm); library(compiler)

# Recoverable
recoverable <- function(simulated)
{

  # Update simulated (if necessary)
  if("original_results" %in% names(simulated)){
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
    cor = "pearson"
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

  # Obtain whether it's recoverable
  if(simulated$parameters$factors == 1){

    # Check if identical
    recoverable <- identical(
      obtain_memberships(aligned$F2, true_membership, simulated$overlap),
      true_membership
    )

  }else{

    # Obtain Omega index
    recoverable <- omega_index(
      obtain_memberships(aligned$F2, true_membership, simulated$overlap),
      true_membership
    ) == 1

  }

  # Return recoverable
  return(list(recoverable = recoverable, max_abs_res = max(abs(EFA$implied_R - simulated$population_correlation))))

}

# Compile function
recoverable <- cmpfun(recoverable)

# Generate condition
generate_condition <- function(condition)
{

  # Total variables
  total_variables <- condition$NVAR * condition$NFAC

  # Set factor sequence
  factor_sequence <- seq_len(condition$NFAC)

  # Number of overlapping cross-loadings to add
  cross_add <- condition$NVAR * 0.25

  # Set unidimensional flag
  unidimensional_flag <- condition$NFAC == 1

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

      # Add overlap
      no_overlap$overlap <- FALSE

      # Check convergence
      convergence <- recoverable(no_overlap)

      # Check if base is recoverable
      if(
        is(no_overlap, "try-error") ||
        !convergence$recoverable ||
        (abs(convergence$max_abs_res) > 0.10)
      ){
        next
      }

      # Check for unidimensional
      if(unidimensional_flag){
        break
      }else{

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

        # Add overlap
        overlap$overlap <- TRUE

        # Check convergence
        convergence <- recoverable(overlap)

        # Check if base is recoverable
        if(
          !is(overlap, "try-error") &&
          convergence$recoverable &&
          (abs(convergence$max_abs_res) < 0.10)
        ){
          break
        }

      }

    }

    # Add error
    no_overlap_error <- add_population_error(
      lf_object = no_overlap, method = "tkl",
      metric = "cfi", fit_value = 0.97
    )

    # Add overlap
    no_overlap_error$overlap <- FALSE

    # Check convergence
    convergence <- recoverable(no_overlap_error)

    # Check if base is recoverable
    if(
      !is(no_overlap_error, "try-error") &&
      convergence$recoverable &&
      (abs(convergence$max_abs_res) < 0.10)
    ){

      # Check for unidimensional
      if(unidimensional_flag){
        break
      }else{

        # Add error
        overlap_error <- add_population_error(
          lf_object = overlap, method = "tkl",
          metric = "cfi", fit_value = 0.97,
          leave_cross_loadings = TRUE
        )

        # Add overlap
        overlap_error$overlap <- TRUE

        # Check convergence
        convergence <- recoverable(overlap_error)

        # Check if base is recoverable
        if(
          !is(overlap_error, "try-error") &&
          convergence$recoverable &&
          (abs(convergence$max_abs_res) < 0.10)
        ){
          break
        }

      }

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

  # Generate samples
  samples <- lapply(c(500, 1000, 2500, 5000), function(N){

    # Generate data
    data <- rmvnorm(n = N, mean = means, sigma = Sigma)

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

}

# Compile function
generate_condition <- cmpfun(generate_condition)
