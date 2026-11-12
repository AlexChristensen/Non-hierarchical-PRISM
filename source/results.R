#%%%%%%%%%%%%%%%%%%%%%%#
#### Obtain Results ####
#%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(psych); library(GPArotation); library(igraph); library(compiler)

# EFA
efa <- function(R, n.obs, nfactors, randomStarts = 100L, ...)
{

  # Data dimensions
  dimensions <- dim(data)

  # Estimate EFA on population correlation matrix
  EFA <- fa(
    r = R, n.obs = n.obs, nfactors = nfactors,
    rotate = "none", fm = "uls"
  )

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

    # Set loadings
    loadings <- rotated$loadings

  }else{

    # Set loadings
    loadings <- EFA$loadings[,seq_len(nfactors)]

  }

  # Return
  return(EGAnet:::reindex_memberships(apply(abs(loadings), 1, which.max)))

}

# Compile
efa <- cmpfun(efa)

# Obtain results
obtain_results <- function(simulated, condition_df)
{

  # Obtain dimensionality
  result <- obtain_methods(simulated)

  # Obtain dimensions
  dimensions <- unlist(result[-c(1:3)])

  # Obtain unique dimensions
  unique_dimensions <- unique(dimensions)

  # Estimate factor memberships
  efa_memberships <- lapply(
    unique_dimensions, function(u){
      efa(
        R = attributes(result)$R,
        n.obs = simulated$parameters$sample_size,
        nfactors = u
      )
    }
  )

  # Set up result
  combined_result <- data.frame(
    METHOD = names(result),
    NDIM = c(
      EGAnet:::unique_length(result$walktrap),
      EGAnet:::unique_length(result$louvain),
      EGAnet:::unique_length(result$prism),
      dimensions
    ),
    ARI = c(
      compare(result$walktrap, simulated$population$membership, method = "adjusted.rand"),
      compare(result$louvain, simulated$population$membership, method = "adjusted.rand"),
      compare(result$prism, simulated$population$membership, method = "adjusted.rand"),
      EGAnet:::nvapply(efa_memberships, compare, simulated$population$membership, method = "adjusted.rand")[
        match(dimensions, unique_dimensions)
      ]
    )
  )

  # Attach condition data frame
  combined_result <- cbind.data.frame(condition_df, combined_result)

  # Return results
  return(combined_result)

}

# Compile function
obtain_results <- cmpfun(obtain_results)