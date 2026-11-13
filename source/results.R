#%%%%%%%%%%%%%%%%%%%%%%#
#### Obtain Results ####
#%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(psych); library(GPArotation); library(igraph); library(compiler)

# EFA
efa <- function(R, n.obs, nfactors, randomStarts = 100L, ...)
{

  # Check for number of factors
  if(nfactors == 0){
    return(NA)
  }

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

# Compile function
efa <- cmpfun(efa)

# Membership comparison
membership_compare <- function(empirical, true)
{

  # Check for missing
  missing <- is.na(empirical)

  # Check for all missing
  if(all(missing)){
    return(NA)
  }

  # Return results
  return(compare(empirical[!missing], true[!missing], method = "adjusted.rand"))

}

# Compile function
membership_compare <- cmpfun(membership_compare)

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
      membership_compare(result$walktrap, simulated$population$membership),
      membership_compare(result$louvain, simulated$population$membership),
      membership_compare(result$prism, simulated$population$membership),
      EGAnet:::nvapply(efa_memberships, membership_compare, simulated$population$membership)[
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