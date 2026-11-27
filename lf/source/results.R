#%%%%%%%%%%%%%%%%%%%%%%#
#### Obtain Results ####
#%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(psych); library(GPArotation); library(igraph); library(compiler)

# Obtain memberships
obtain_memberships <- function(loadings, true_membership, overlap)
{

  # Check for NA
  if(all(is.na(loadings))){
    return(NA)
  }

  # Obtain memberships
  factor_sequence <- seq_len(ncol(loadings))
  primary <- outer(max.col(loadings), factor_sequence, "==")

  # Check for overlap
  if(overlap){
    secondary <- outer(max.col(abs(loadings * (!primary))), factor_sequence, "==")
    secondary[rowSums(true_membership) == 1,] <- FALSE
    return(primary | secondary)
  }else{
    return(primary)
  }

}

# Compile
obtain_memberships <- cmpfun(obtain_memberships)

# Obtain omega index
# Collins, L. M., & Dent, C. W. (1988).
# Omega: A general formulation  of the rand index of cluster recovery suitable for non-disjoint solutions.
# Multivariate Behavioral Research, 23(2), 231–242.
omega_index <- function(membership, true_membership)
{

  # Check for NA
  if(all(is.na(membership))){
    return(NA)
  }

  # Co-membership counts
  comembership <- tcrossprod(membership)
  true_comembership <- tcrossprod(true_membership)

  # Obtain lower triangle
  lower_triangle <- lower.tri(comembership)
  lower_comembership <- comembership[lower_triangle]
  true_lower_comembership <- true_comembership[lower_triangle]

  # Guard against single dimensions
  if(length(unique(true_lower_comembership)) == 1){
    return(NA)
  }

  # Total comparisons
  total <- length(lower_comembership)

  # Maximum overlapping communities (used for bins below)
  max_overlap <- max(lower_comembership, true_lower_comembership, na.rm = TRUE) + 1

  # omega_u: fraction of pairs whose co-occurrence count matches exactly
  omega_u <- mean(lower_comembership == true_lower_comembership, na.rm = TRUE)

  # omega_e: expected match fraction under independent marginal distributions
  freq <- tabulate(lower_comembership + 1, nbins = max_overlap) / total
  true_freq <- tabulate(true_lower_comembership + 1, nbins = max_overlap) / total
  omega_e <- sum(freq * true_freq)

  # Return omega index
  return((omega_u - omega_e) / (1 - omega_e))

}

# Compile
omega_index <- cmpfun(omega_index)

# EFA
efa <- function(data, nfactors, corr, n.obs = NULL, randomStarts = 100L, ...)
{

  # Data dimensions
  dimensions <- dim(data)

  # Check for correlations
  if(dimensions[1] != dimensions[2]){
    data <- switch(
      corr,
      "pearson" = cor(data),
      "polychoric" = polychoric(data)
    )
    n.obs <- dimensions[1]
  }else if(is.null(n.obs)){
    n.obs <- 8.3e09
  }

  # Estimate EFA on population correlation matrix
  EFA <- fa(
    r = data, n.obs = n.obs, nfactors = nfactors,
    rotate = "none", fm = EGAnet:::swiftelse(corr == "pearson", "ml", "uls")
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

    # Update output
    EFA$loadings <- rotated$loadings
    EFA$Phi <- rotated$Phi

    # Obtain model-implied correlations
    EFA$implied_R <- EFA$loadings %*% tcrossprod(EFA$Phi, EFA$loadings)
    diag(EFA$implied_R) <- 1

  }else{

    # Obtain model-implied correlations
    EFA$Phi <- matrix(1)
    EFA$implied_R <- tcrossprod(EFA$loadings)
    diag(EFA$implied_R) <- 1

  }

  # Return
  return(EFA)

}

# Compile
efa <- cmpfun(efa)

# EFA loadings
efa_loadings <- function(R, n.obs, nfactors)
{

  # Check for number of factors
  if(is.na(nfactors) || nfactors < 1){
    return(NA)
  }

  # Estimate EFA (guard against too many factors)
  EFA <- try(
    EGAnet:::silent_call(efa(
      data = R, nfactors = nfactors, corr = "polychoric", n.obs = n.obs
    )), silent = TRUE
  )

  # Check for error
  if(is(EFA, "try-error")){
    return(NA)
  }

  # Return memberships
  return(abs(EFA$loadings))

}

# Compile function
efa_loadings <- cmpfun(efa_loadings)

# Obtain results
obtain_results <- function(samples, condition_df)
{

  # Set error flag
  error_flag <- "original_results" %in% names(samples)

  # Set up true membership
  if(error_flag){
    true_membership <- abs(samples$original_results$parameters$loadings) >= 0.40
  }else{
    true_membership <- abs(samples$parameters$loadings) >= 0.40
  }

  # Obtain dataset
  return(
    do.call(rbind.data.frame, lapply(samples$data, function(data){

      # Obtain sample size
      N <- dim(data)[1]

      # Obtain dimensionality
      result <- EGAnet:::silent_call(obtain_methods(data))

      # Obtain dimensions
      dimensions <- unlist(result[-c(1:3)])

      # Obtain unique dimensions
      unique_dimensions <- unique(dimensions)

      # Obtain network membership
      ega_memberships <- lapply(
        result[1:3], obtain_memberships,
        true_membership = true_membership,
        overlap = samples$overlap
      )

      # Estimate factor memberships
      efa_memberships <- lapply(
        unique_dimensions, function(u){
          obtain_memberships(
            loadings = efa_loadings(
              R = attributes(result)$R,
              n.obs = N,
              nfactors = u
            ), true_membership = true_membership,
            overlap = samples$overlap
          )
        }
      )

      # Set up result
      combined_result <- data.frame(
        METHOD = names(result),
        N = N, OVERLAP = samples$overlap, ERROR = error_flag,
        NDIM = c(
          EGAnet:::nvapply(result[c(1:3)], ncol),
          dimensions
        ),
        ARI = c(
          EGAnet:::nvapply(ega_memberships, omega_index, true_membership),
          EGAnet:::nvapply(efa_memberships, omega_index, true_membership)[
            match(dimensions, unique_dimensions)
          ]
        )
      )

      # Attach condition data frame
      return(cbind.data.frame(condition_df, combined_result))

    }))
  )

}

# Compile function
obtain_results <- cmpfun(obtain_results)
