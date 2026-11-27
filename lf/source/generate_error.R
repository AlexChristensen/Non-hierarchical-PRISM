#%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Generate Condition ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(latentFactoR); library(fungible); library(psych); library(compiler)

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
  EFA <- try(
    EGAnet:::silent_call(efa(
      data = simulated$population_correlation,
      nfactors = simulated$parameters$factors,
      corr = "pearson"
    )), silent = TRUE
  )

  # Check for estimation failure
  if(is(EFA, "try-error")){
    return(list(recoverable = FALSE, max_abs_res = NA))
  }

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

  # Calculate implied R
  implied_R <- EFA$loadings %*% tcrossprod(EFA$Phi, EFA$loadings)
  diag(implied_R) <- 1

  # Set true membership
  true_membership <- abs(true_loadings) >= 0.40

  # Obtain whether it's recoverable
  recoverable <- omega_index(
    obtain_memberships(aligned$F2, true_membership, simulated$overlap),
    true_membership
  ) == 1

  # Return recoverable
  return(list(recoverable = recoverable, max_abs_res = max(abs(implied_R - simulated$population_correlation))))

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

      # Check if base is recoverable
      if(is(no_overlap, "try-error") || !recoverable(no_overlap)$recoverable){
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

      # Add overlap
      overlap$overlap <- TRUE

      # Check if base is recoverable
      if(!is(overlap, "try-error") && recoverable(overlap)$recoverable){
        break
      }

    }

    # Add error to no overlap
    no_overlap_result <- do.call(
      rbind.data.frame, lapply(c("cb", "tkl", "wb", "hb"), function(method){

        # Switch metric based on method
        metric <- switch(
          method,
          "cb" = "cfi",
          "tkl" = "cfi",
          "wb" = "rmsea",
          "hb" = "srmr"
        )

        # Switch value based on metric
        values <- switch(
          metric,
          "cfi" = seq(0.99, 0.95, -0.01),
          "rmsea" = seq(0.01, 0.05, 0.01),
          "srmr" = seq(0.01, 0.05, 0.01)
        )

        # Loop over values
        fit_values <- do.call(
          rbind, lapply(values, function(value){

            # Add population error
            no_overlap_error <- latentFactoR::add_population_error(
              lf_object = no_overlap, method = method,
              metric = metric, fit_value = value
            )

            # Add overlap
            no_overlap_error$overlap <- FALSE

            # Return fit measures
            return(
              c(
                method = method,
                metric = metric,
                value = value,
                no_overlap_error$population_error$fit,
                max_residual = no_overlap_error$population_error$max_error,
                recoverable = recoverable(no_overlap_error)$recoverable
              )
            )

          })
        )
      })
    )

    # Add error to overlap
    overlap_result <- do.call(
      rbind.data.frame, lapply(c("cb", "tkl", "wb", "hb"), function(method){

        # Switch metric based on method
        metric <- switch(
          method,
          "cb" = "cfi",
          "tkl" = "cfi",
          "wb" = "rmsea",
          "hb" = "srmr"
        )

        # Switch value based on metric
        values <- switch(
          metric,
          "cfi" = seq(0.99, 0.95, -0.01),
          "rmsea" = seq(0.01, 0.05, 0.01),
          "srmr" = seq(0.01, 0.05, 0.01)
        )

        # Loop over values
        fit_values <- do.call(
          rbind, lapply(values, function(value){

            # Add population error
            overlap_error <- latentFactoR::add_population_error(
              lf_object = overlap, method = method,
              metric = metric, fit_value = value,
              leave_cross_loadings = TRUE
            )

            # Add overlap
            overlap_error$overlap <- TRUE

            # Return fit measures
            return(
              c(
                method = method,
                metric = metric,
                value = value,
                overlap_error$population_error$fit,
                max_residual = overlap_error$population_error$max_error,
                recoverable = recoverable(overlap_error)$recoverable
              )
            )

          })
        )
      })
    )

    # All clear
    break

  }

  # Construct data frame
  no_overlap <- as.data.frame(no_overlap_result)
  no_overlap$overlap <- FALSE
  overlap <- as.data.frame(overlap_result)
  overlap$overlap <- TRUE

  # Return data
  return(rbind.data.frame(no_overlap, overlap))

}

# Compile function
generate_condition <- cmpfun(generate_condition)
