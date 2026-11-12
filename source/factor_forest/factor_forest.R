#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Factor Forest (Goretzko) ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Estimates the number of dimensions (factors) with the pre-trained xgboost
# model of Goretzko & Buhner (2020; https://osf.io/mvrau/), applied to ordinal
# data as in Goretzko & Buhner (2022; Applied Psychological Measurement).
#
# Feature extraction follows `analyze()` in `Evaluation.R` (the authors' own
# evaluation code for ordinal data). The model's 184 features are:
#   N, p, 22 summary features of the correlation matrix and its eigenvalues,
#   parallel analysis, empirical Kaiser criterion, comparison data approach,
#   and the first 80 eigenvalues of the correlation matrix and of the
#   reduced correlation matrix (FA).
#
# For ordinal data (the default), polychoric correlations replace Pearson
# correlations for the features, parallel analysis, and EKC. The comparison
# data approach runs on the raw data with Pearson correlations (as in
# `Evaluation.R`).
#
# The model is `tunedxgb.json`, a one-time conversion of `tunedxgb.rds` that
# current xgboost versions can read (see `convert_model.md`).
#
# Usage (from the project root, as in the other simulation scripts):
#   source("./source/factor_forest/factor_forest.R")
#   factor_forest(simulated$data)
#   factor_forest(simulated$data, R = R) # re-use polychoric matrix already computed

# Load packages
library(xgboost); library(psych); library(ineq); library(compiler)

# Locate the folder holding this script, the CD approach, and the model
.factor_forest_dir <- local({

  # Candidate folders (explicit option > location of this file > convention)
  candidates <- c(
    getOption("factor_forest.dir"),
    tryCatch(dirname(sys.frame(1)$ofile), error = function(e){NULL}),
    "./source/factor_forest", "./factor_forest", "."
  )

  # Return the first with the model in it
  for(path in candidates){
    if(file.exists(file.path(path, "tunedxgb.json"))){
      return(path)
    }
  }

  stop(
    "Could not find `tunedxgb.json`. Set `options(factor_forest.dir = \"<path>\")`",
    call. = FALSE
  )

})

# Dimensionality methods used for the features (`source/methods.R`):
# `kaiser_guttman()`, `empirical_kaiser()`, `parallel_analysis()`, and
# `comparison_data()`; loaded here unless already available
if(!all(vapply(
  c("kaiser_guttman", "empirical_kaiser", "parallel_analysis", "comparison_data"),
  exists, logical(1)
))){
  source(Filter(
    file.exists,
    c(
      file.path(.factor_forest_dir, "..", "methods.R"), "./source/methods.R"
    )
  )[1])
}

# Cache so the model is only read once per R process
.factor_forest_cache <- new.env(parent = emptyenv())

# Load (and cache) the pre-trained model
factor_forest_model <- function(path = file.path(.factor_forest_dir, "tunedxgb.json"))
{

  # Read once
  if(is.null(.factor_forest_cache[[path]])){
    .factor_forest_cache[[path]] <- xgboost::xgb.load(path)
  }

  # Return model
  return(.factor_forest_cache[[path]])

}

# Number of eigenvalues the model uses as features
.factor_forest_eigenvalues <- 80L

# Pad (with -1000, as in training) or truncate to the model's expected length
.factor_forest_pad <- function(x, length = .factor_forest_eigenvalues)
{

  # Initialize padded
  padded <- rep(-1000, length)

  # Fill
  n <- min(base::length(x), length)
  padded[seq_len(n)] <- x[seq_len(n)]

  # Return padded
  return(padded)

}

# Compute the model's features (one-row data frame)
#
# `data`: N x p matrix or data frame of item responses
# `R`: optional p x p correlation matrix (e.g., polychoric already computed);
#      computed from `data` when NULL
# `cor`: "poly" (ordinal data; default) or "pearson" (continuous data)
factor_forest_features <- function(
    data, R = NULL, cor = c("poly", "pearson"),
    use = "pairwise.complete.obs", F.Max = 8
)
{

  # Check for correlation type
  cor <- match.arg(cor)

  # Ensure numeric matrix
  data <- as.matrix(data)

  # Obtain dimensions
  N <- nrow(data)
  p <- ncol(data)

  # Correlation matrix
  if(is.null(R)){
    dat_cor <- switch(
      cor,
      "poly" = psych::polychoric(data)$rho,
      "pearson" = stats::cor(data, use = use)
    )
  }else{
    dat_cor <- as.matrix(R)
  }

  # Eigenvalues
  eigval <- eigen(dat_cor, symmetric = TRUE)$values
  vareig <- cumsum(eigval) / p

  # Empirical Kaiser criterion
  ekc <- empirical_kaiser(N = N, eigenvalues = eigval)

  # Eigenvalue based features (including the Kaiser-Guttman rule)
  eiggreater1 <- kaiser_guttman(eigenvalues = eigval)
  releig1 <- eigval[1] / p
  releig2 <- sum(eigval[1:2]) / p
  releig3 <- sum(eigval[1:3]) / p
  eiggreater07 <- sum(eigval > 0.7)
  sdeigval <- stats::sd(eigval)
  var50 <- min(which(vareig > 0.50))
  var75 <- min(which(vareig > 0.75))

  # Matrix norm features
  onenorm <- norm(dat_cor, "O")
  frobnorm <- norm(dat_cor, "F")
  maxnorm <- norm(dat_cor - diag(p), "M")
  avgcor <- sum(abs(dat_cor - diag(p))) / (p * (p - 1))
  specnorm <- sqrt(eigen(crossprod(dat_cor), symmetric = TRUE)$values[1])

  smlcor <- sum(dat_cor <= 0.1)
  avgcom <- mean(suppressWarnings(psych::smc(dat_cor)))
  det <- det(dat_cor)

  # "Inequality" features
  #
  # NOTE: the model was trained with these computed on `lower.tri(dat_cor)`
  # (a logical matrix) rather than on the correlations themselves (see
  # `Evaluation.R`). This is kept on purpose: features must be computed exactly
  # as they were in training for the pre-trained model's predictions to be valid.
  KMO <- suppressWarnings(psych::KMO(dat_cor)$MSA)
  Gini <- ineq::ineq(lower.tri(dat_cor), type = "Gini")
  Kolm <- ineq::ineq(lower.tri(dat_cor), type = "Kolm")

  # Parallel analysis (common factor model, as in training)
  pa <- parallel_analysis(
    data, corr = switch(cor, "poly" = "polychoric", "pearson" = "pearson"),
    R = dat_cor, type = "fa"
  )
  pa_solution <- pa[["fa"]]
  fa_eigval <- attr(pa, "eigenvalues")$fa

  # Pad (or truncate) eigenvalues to the model's expected length
  eigval <- .factor_forest_pad(eigval)
  fa_eigval <- .factor_forest_pad(fa_eigval)
  names(eigval) <- paste0("eigval", seq_along(eigval))
  names(fa_eigval) <- paste0("fa_eigval", seq_along(fa_eigval))

  # Comparison data approach
  cd <- comparison_data(Data = data, F.Max = F.Max, use = use)

  # Return features (column order is the order the model was trained on)
  return(
    cbind(
      data.frame(
        N = N, p = p, eiggreater1, releig1, releig2, releig3, eiggreater07,
        sdeigval, var50, var75, onenorm, frobnorm, maxnorm, avgcor, specnorm,
        smlcor, avgcom, det, KMO, Gini, Kolm, pa_solution, ekc, cd
      ),
      t(eigval), t(fa_eigval)
    )
  )

}

# Compile function
factor_forest_features <- cmpfun(factor_forest_features)

# Estimate the number of factors
#
# Returns an integer between 1 and 8 (the model's class labels)
factor_forest <- function(
    data, R = NULL, model = NULL, cor = c("poly", "pearson"),
    use = "pairwise.complete.obs", F.Max = 8
)
{

  # Obtain (cached) pre-trained model
  if(is.null(model)){
    model <- factor_forest_model()
  }

  # Compute features
  features <- factor_forest_features(
    data, R = R, cor = cor, use = use, F.Max = F.Max
  )

  # Predict class probabilities (uses the best early-stopping iteration)
  probabilities <- matrix(
    suppressWarnings(predict(model, as.matrix(features))),
    nrow = 1
  )

  # Return number of factors (class labels are 1, ..., 8)
  return(max.col(probabilities, ties.method = "first"))

}

# Compile function
factor_forest <- cmpfun(factor_forest)
