#%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Generate Condition ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(psychNetwoRk); library(EGAnet); library(igraph)
library(mvtnorm); library(compiler)

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
  R <- simulated$population$R; diag(R) <- 0
  unidimensional <- community.consensus(R)
  walktrap <- community.detection(Omega)
  louvain <- community.consensus(Omega)
  prism <- community.prism(Omega)

  # Return recoverable
  return(
    EGAnet:::swiftelse(
      simulated$parameters$blocks == 1,
      EGAnet:::unique_length(unidimensional) == 1,
      any(
        c(
          walktrap = ifelse(
            any(is.na(walktrap)), FALSE,
            abs(compare(walktrap, simulated$population$membership, method = "adjusted.rand") - 1) < 1e-06
          ),
          louvain = ifelse(
            any(is.na(louvain)), FALSE,
            abs(compare(louvain, simulated$population$membership, method = "adjusted.rand") - 1) < 1e-06
          ),
          prism = ifelse(
            any(is.na(prism)), FALSE,
            abs(compare(prism, simulated$population$membership, method = "adjusted.rand") - 1) < 1e-06
          )
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
  max_residual <- simulated$population_error$max_error < 0.10
  max_discrepancy <- max(abs(simulated$population$discrepancy)) < 0.10

  # Unidimensional: CFI is limited by the few degrees of freedom
  # so only the residual and discrepancy criteria are used
  if(simulated$parameters$blocks == 1){
    return(all(max_residual, max_discrepancy))
  }

  # Multidimensional
  cfi <- abs(simulated$population_error$fit[["CFI"]] - 0.97) < 0.01

  # Return convergence
  return(all(cfi, max_residual, max_discrepancy))

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

    # Generate until recoverable
    while(TRUE){

      # Generate data
      simulated_base <- try(simulate_sbm(
        nodes = condition$NODES, blocks = condition$BLOCKS,
        density_matrix = DM, diffusion = condition$DIFFUSION,
        sample_size = 1, negative_proportion = 0
      ), silent = TRUE)

      # Check if base is recoverable
      if(!is(simulated_base, "try-error") && recoverable(simulated_base)){
        break
      }

    }

    # Set CFI target
    target <- 0.97

    # Set CFI target
    if(condition$BLOCKS == 1){

      # Loop by increasing CFI target
      while(TRUE){

        # Obtain error data
        simulated_error <- try(add_population_error(
          ggm_object = simulated_base, method = "tkl", metric = "cfi",
          fit_value = target
        ), silent = TRUE)

        # Check if base is recoverable
        if(is(simulated_error, "try-error") || convergence(simulated_error)){
          break
        }

        # Increase target
        target <- sqrt(target)

      }

    }

    # Obtain error data
    simulated_error <- try(add_population_error(
      ggm_object = simulated_base, method = "tkl", metric = "cfi",
      fit_value = target
    ), silent = TRUE)

    # Check if base is recoverable
    if(!is(simulated_error, "try-error") && convergence(simulated_error) && recoverable(simulated_error)){
      break
    }

  }

  # Obtain number of variables
  n_variables <- length(simulated_base$population$membership)

  # Set skew
  skewness <- psychNetwoRk:::runif_xoshiro(n_variables, min = 0, max = 2)

  # Obtain Cholesky decomposition
  cholesky_base <- chol(simulated_base$population$R)
  cholesky_error <- chol(simulated_error$population$R)

  # Obtain means and Sigma
  means <- rep(0, n_variables); Sigma <- diag(n_variables)

  # Generate samples
  samples <- lapply(c(500, 1000, 2500, 5000), function(N){

    # Generate data
    data <- rmvnorm(n = N, mean = means, sigma = Sigma)

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


  # Return data
  return(samples)

}

# Compile function
generate_condition <- cmpfun(generate_condition)
