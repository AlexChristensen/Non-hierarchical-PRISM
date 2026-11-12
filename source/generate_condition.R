#%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Generate Condition ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(psychNetwoRk); library(EGAnet); library(igraph); library(compiler)

# Recoverable
recoverable <- function(simulated)
{

  # Check for pseudo Omega
  if(is.null(simulated$population$Omega_pseudo)){
    Omega <- simulated$population$Omega
  }else{
    Omega <- simulated$population$Omega_pseudo
  }

  # Estimate memberships
  walktrap <- community.detection(Omega)
  louvain <- community.consensus(Omega)
  prism <- community.prism(Omega)

  # Return recoverable
  return(
    any(
      c(
        walktrap = ifelse(
          any(is.na(walktrap)), FALSE,
          abs(compare(walktrap, simulated$population$membership, method = "adjusted.rand") - 1) < 1e06
        ),
        louvain = ifelse(
          any(is.na(louvain)), FALSE,
          abs(compare(louvain, simulated$population$membership, method = "adjusted.rand") - 1) < 1e06
        ),
        prism = ifelse(
          any(is.na(prism)), FALSE,
          abs(compare(prism, simulated$population$membership, method = "adjusted.rand") - 1) < 1e06
        )
      )
    )
  )

}

# Compile function
recoverable <- cmpfun(recoverable)

# Convergence Criterion
convergence <- function(simulated)
{

  # Set up logicals
  cfi <- simulated$population_error$fit[["CFI"]] > 0.95
  rmsea <- abs(simulated$population_error$fit[["RMSEA"]] - 0.02) < 1
  max_residual <- simulated$population_error$max_error < 0.10

  # Return convergence
  return(all(cfi, rmsea, max_residual))

}

# Compile function
convergence <- cmpfun(convergence)

# Generate condition
generate_condition <- function(condition)
{

  # Set density matrix
  DM <- matrix(condition$BETWEEN, nrow = condition$BLOCKS, ncol = condition$BLOCKS)
  diag(DM) <- condition$WITHIN

  # Generate until recoverable
  while(TRUE){

    # Generate data
    simulated_base <- psychNetwoRk:::silent_call(simulate_sbm(
      nodes = condition$NODES, blocks = condition$BLOCKS,
      density_matrix = DM, diffusion = condition$DIFFUSION,
      sample_size = condition$N
    ))

    # Check if base is recoverable
    if(recoverable(simulated_base)){
      break
    }

  }

  # Generate until recoverable
  while(TRUE){

    # Obtain error data
    simulated_error <- psychNetwoRk:::silent_call(add_population_error(
      ggm_object = simulated_base, method = "wb", metric = "rmsea",
      fit_value = switch(
        as.character(simulated_base$parameters$blocks),
        "2" = 0.04, "4" = 0.03, "8" = 0.02
      )
    ))

    # Check if base is recoverable
    if(convergence(simulated_error) && recoverable(simulated_error)){
      break
    }

  }

  # Obtain number of variables
  n_variables <- length(simulated_base$population$membership)

  # Set skew
  skewness <- psychNetwoRk:::runif_xoshiro(n_variables, min = 0, max = 2)

  # Add skew to condition
  simulated_error$parameters$skew <- simulated_base$parameters$skew <- skewness

  # Add skew to the data
  for(v in seq_len(n_variables)){

    # Add skew
    ## Base
    simulated_base$data[,v] <- categorize(
      simulated_base$data[,v], categories = 5, skew_value = skewness[v]
    )

    ## Error
    simulated_error$data[,v] <- categorize(
      simulated_error$data[,v], categories = 5, skew_value = skewness[v]
    )

  }


  # Return data
  return(list(no_error = simulated_base, error = simulated_error))

}

# Compile function
generate_condition <- cmpfun(generate_condition)
