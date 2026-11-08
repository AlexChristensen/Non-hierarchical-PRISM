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

}

# Compile function
recoverable <- cmpfun(recoverable)

# Generate condition
generate_condition <- function(condition)
{

  # Set density matrix
  DM <- matrix(condition$BETWEEN, nrow = condition$BLOCKS, ncol = condition$BLOCKS)
  diag(DM) <- condition$WITHIN

  # Generate until recoverable
  while(TRUE){

    # Generate data
    simulated_base <- simulate_sbm(
      nodes = condition$NODES, blocks = condition$BLOCKS,
      density_matrix = DM, diffusion = condition$DIFFUSION,
      sample_size = condition$N
    )

    # Check if base is recoverable
    if(any(recoverable(simulated_base))){
      break
    }

  }

  # Loop over error method
  return(
    do.call(
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
            simulated_error <- add_population_error(
              ggm_object = simulated_base, method = method,
              metric = metric, fit_value = value
            )

            # Return fit measures
            return(
              c(
                method = method,
                metric = metric,
                value = value,
                simulated_error$population_error$fit,
                max_residual = simulated_error$population_error$max_error,
                recoverable(simulated_error)
              )
            )

          })
        )

      })
    )
  )

}

# Compile function
generate_condition <- cmpfun(generate_condition)
