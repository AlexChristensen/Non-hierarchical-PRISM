#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Dimensionality Methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(EGAnet); library(latentFactoR); library(psych)
library(GPArotation); library(compiler); library(Matrix)

# Set polychoric
polychoric <- function(data, empty.method = "zero", empty.value = "point_five")
{

  # Obtain polychoric matrix
  R <- polychoric.matrix(data, empty.method = empty.method, empty.value = empty.value)

  # Check for smoothing
  if(any(eigen(R, symmetric = TRUE, only.values = TRUE)$values < 0)){
    R <- as.matrix(nearPD(R, corr = TRUE, keepDiag = TRUE, ensureSymmetry = TRUE)$mat)
  }

  # Return correlation matrix
  return(R)

}

# Compile function
polychoric <- cmpfun(polychoric)

# Parallel Analysis
#
# `type` selects which variants are computed (default: all three):
#   "pca"    = principal components vs. mean eigenvalues of shuffled data
#   "paf"    = principal axis factoring (SMC communalities) vs. mean
#              eigenvalues of shuffled data
#   "paf_fa" = common factor model (minres) vs. the `quant` quantile of
#              eigenvalues from random normal data, as in
#              `psych::fa.parallel(fa = "fa")`. This is the variant used by
#              `factor_forest`. Its reference only depends on N and p, so
#              (unlike `psych::fa.parallel`) no correlations of resampled data
#              are computed, which makes it much faster. Settings match the
#              original `fa.parallel(fa = "fa", cor = "poly")` call: random
#              normal data (psych's default `sim = TRUE`; shuffled variables
#              are only used for `nfact` when `sim = FALSE`), 20 iterations,
#              95th percentile, one-factor minres eigenvalues
# `R` is an optional, already computed empirical correlation matrix
# The empirical eigenvalues of each requested type are attached as the
# attribute "eigenvalues"
parallel_analysis <- function(
    data, corr, n.iter = 20, quant = 0.95, R = NULL, type = c("pca", "paf", "paf_fa")
)
{

  # Check for type
  type <- match.arg(type, choices = c("pca", "paf", "paf_fa"), several.ok = TRUE)

  # Set up correlation function
  corr_FUN <- switch(
    corr,
    "pearson" = cor,
    "polychoric" = polychoric
  )

  # Sub-routine for correct counting of dimensions
  count_dimensions <- function(greater_than, variables)
  {

    # Initialize count
    count <- 0

    # Loop over to count consecutive dimensions
    for(i in seq_len(variables)){

      # Check for element TRUE
      if(greater_than[i]){
        count <- count + 1
      }else{
        break
      }

    }

    # Return dimensions
    return(count)

  }

  # Set dimensions
  dimensions <- dim(data)
  variables <- dimensions[2]

  # Compute (or use supplied) empirical correlations
  if(is.null(R)){
    R_empirical <- corr_FUN(data)
  }else{
    R_empirical <- as.matrix(R)
  }

  # Keep a copy (the diagonal is replaced below for principal axis factoring)
  R_correlation <- R_empirical

  # Initialize results
  counts <- list()
  eigenvalues <- list()

  # Shuffled data variants
  if(any(c("pca", "paf") %in% type)){

    # Obtain PCA empirical eigenvalues
    PCA_empirical <- eigen(R_empirical, symmetric = TRUE, only.values = TRUE)$values

    # Compute SMC
    diag(R_empirical) <- psych::smc(R_empirical)

    # Obtain PAF empirical eigenvalues
    PAF_empirical <- eigen(R_empirical, symmetric = TRUE, only.values = TRUE)$values

    # Collect empirical resamples
    resampled <- lapply(seq_len(n.iter), function(i){

      # Shuffle data
      shuffled <- apply(data, 2, psychNetwoRk:::shuffle)

      # Compute polychoric correlations
      R_shuffled <- corr_FUN(shuffled)

      # Obtain PCA shuffled eigenvalues
      PCA_shuffled <- eigen(R_shuffled, symmetric = TRUE, only.values = TRUE)$values

      # Compute SMC
      diag(R_shuffled) <- psych::smc(R_shuffled)

      # Return eigenvalues
      return(list(
        pca = PCA_shuffled,
        paf = eigen(R_shuffled, symmetric = TRUE, only.values = TRUE)$values
      ))

    })

    # Obtain means
    PCA_mean <- rowMeans(do.call(cbind, lapply(resampled, function(x){x$pca})), na.rm = TRUE)
    PAF_mean <- rowMeans(do.call(cbind, lapply(resampled, function(x){x$paf})), na.rm = TRUE)

    # # Obtain quantiles
    # PCA_quantile <- apply(
    #   do.call(cbind, lapply(resampled, function(x){x$pca})), 1, quantile, probs = quant
    # )
    # PAF_quantile <- apply(
    #   do.call(cbind, lapply(resampled, function(x){x$paf})), 1, quantile, probs = quant
    # )

    # Store results
    counts$pca <- count_dimensions(PCA_empirical >= PCA_mean, variables)
    counts$paf <- count_dimensions(PAF_empirical >= PAF_mean, variables)
    eigenvalues$pca <- PCA_empirical
    eigenvalues$paf <- PAF_empirical

  }

  # Common factor model (minres) vs. random normal data (as `psych::fa.parallel`)
  if("paf_fa" %in% type){

    # Obtain empirical common factor eigenvalues
    FA_empirical <- suppressWarnings(suppressMessages(
      psych::fa(
        R_correlation, nfactors = 1, rotate = "none",
        fm = "minres", warnings = FALSE
      )$values
    ))

    # Obtain common factor eigenvalues of random normal data
    FA_random <- vapply(seq_len(n.iter), function(i){

      # Compute correlations of random normal data
      R_random <- cor(
        matrix(psychNetwoRk:::rnorm_ziggurat(dimensions[1] * variables), nrow = dimensions[1])
      )

      # Return eigenvalues
      return(
        suppressWarnings(suppressMessages(
          psych::fa(
            R_random, nfactors = 1, rotate = "none",
            fm = "minres", SMC = FALSE, warnings = FALSE
          )$values
        ))
      )

    }, numeric(variables))

    # Obtain quantiles
    FA_quantile <- apply(FA_random, 1, quantile, probs = quant)

    # Store results
    counts$paf_fa <- count_dimensions(FA_empirical > FA_quantile, variables)
    eigenvalues$paf_fa <- FA_empirical

  }

  # Collect results (in the order requested)
  results <- unlist(counts[type])
  attr(results, "eigenvalues") <- eigenvalues[type]

  # Return results
  return(results)

}

# Compile
parallel_analysis <- cmpfun(parallel_analysis)

# Kaiser-Guttman rule: number of eigenvalues greater than one
#
# `R` is a correlation matrix; supply `eigenvalues` to skip the decomposition
kaiser_guttman <- function(R = NULL, eigenvalues = NULL)
{

  # Obtain eigenvalues
  if(is.null(eigenvalues)){
    eigenvalues <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  }

  # Return dimensions
  return(sum(eigenvalues > 1))

}

# Compile
kaiser_guttman <- cmpfun(kaiser_guttman)

# Empirical Kaiser criterion (Braeken & van Assen, 2017)
#
# Restricted version: reference eigenvalues are never below one
# `R` is a correlation matrix and `N` is the sample size;
# supply `eigenvalues` to skip the decomposition
empirical_kaiser <- function(R = NULL, N, eigenvalues = NULL)
{

  # Obtain eigenvalues
  if(is.null(eigenvalues)){
    eigenvalues <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  }

  # Set number of variables
  variables <- length(eigenvalues)

  # Compute reference eigenvalues
  reference <- rep(0, variables)
  for(i in seq_len(variables)){
    reference[i] <- max(
      ((1 + sqrt(variables / N))^2) * (variables - sum(reference)) / (variables - i + 1),
      1
    )
  }

  # Return dimensions (eigenvalues before the first one not exceeding reference)
  return(which(eigenvalues <= reference)[1] - 1)

}

# Compile
empirical_kaiser <- cmpfun(empirical_kaiser)

#### Comparison Data ####
#
# Optimized implementation of the comparison data approach of Ruscio and Roche
# (2012). It is a drop-in replacement for `EFA.Comp.Data()` in
# `factor_forest/CD_approach.R` (the original, kept as the reference): same
# algorithm, defaults, and random number stream, but much faster.
#
# What was changed (no change to the algorithm):
#   - Correlations of resampled data use BLAS (`crossprod`) rather than
#     `cor(use = "pairwise.complete.obs")` (~8x faster for the resampling
#     loop). The comparison population is built by resampling observed
#     non-missing values, so it never has missing values.
#   - Only eigenvalues are computed where eigenvectors are not needed, and
#     `eigen(symmetric = TRUE)` skips the symmetry check.
#   - The shared + unique component matrix is computed once per iteration
#     rather than once per variable.
#   - Imposing the observed distributions no longer physically re-sorts the
#     whole population matrix once per variable; the same row permutation is
#     tracked as an index vector and applied once.
#   - Values that do not change across iterations (per-variable distributions,
#     target correlation matrix and its eigenvalues) are computed once.
#
# Checked against `EFA.Comp.Data()`: with the original's `cor()` call the
# generated populations are identical. The BLAS correlation differs from
# `cor()` at the ~1e-17 level, which the search over intermediate correlation
# matrices can amplify when eigenvalues are nearly tied (e.g., symmetric
# blocks), so the final count can differ from the original's for the same seed
# in a minority of data sets (the original differs from itself across seeds
# about as often). The `Graph` argument was dropped.

# Correlation matrix of data that has no missing values
.cd_cor <- function(x, cor_type)
{

  # Spearman needs the ranking done by `cor()`
  if(cor_type == "spearman"){
    return(stats::cor(x, method = "spearman"))
  }

  # Pearson via BLAS
  x <- x - rep(colMeans(x), each = nrow(x))
  cross <- crossprod(x)
  scale <- sqrt(diag(cross))

  # Return correlations
  return(cross / tcrossprod(scale))

}

# Principal axis factor analysis on a correlation matrix (loadings only)
.cd_factor_analysis <- function(cor_matrix, n_factors, max_iter = 50, criterion = 0.001)
{

  # Initialize
  k <- ncol(cor_matrix)
  old_h2 <- rep(99, k)
  change <- 1
  iter <- 0

  # Iterate communalities
  while((change >= criterion) && (iter < max_iter)){

    # Next iteration
    iter <- iter + 1

    # Eigendecomposition
    eig <- eigen(cor_matrix, symmetric = TRUE)

    # Loadings and communalities
    loadings <- eig$vectors[, seq_len(n_factors), drop = FALSE] *
      rep(sqrt(eig$values[seq_len(n_factors)]), each = k)
    h2 <- rowSums(loadings * loadings)

    # Update
    change <- max(abs(old_h2 - h2))
    old_h2 <- h2
    diag(cor_matrix) <- h2

  }

  # Return loadings
  return(loadings)

}

# Generate comparison population
.cd_generate <- function(
    data, n_factors, N, cor_type, target_cor, distribution_values,
    max_trials = 5, initial_multiplier = 1
)
{

  # Initialize
  k <- ncol(data)
  best_rmsr <- 1
  trials_without_improvement <- 0
  lower <- lower.tri(target_cor)
  n_lower <- sum(lower)

  # Distribution for each variable (sorted bootstrap of observed values)
  distributions <- matrix(0, nrow = N, ncol = k)
  for(i in seq_len(k)){
    distributions[, i] <- sort(
      sample(distribution_values[[i]], size = N, replace = TRUE)
    )
  }

  # Shared and unique components
  shared_comp <- matrix(rnorm(N * n_factors, 0, 1), nrow = N, ncol = n_factors)
  unique_comp <- matrix(rnorm(N * k, 0, 1), nrow = N, ncol = k)

  # Loadings -> data with the observed distributions imposed
  #
  # Rows are returned in the same order the original implementation produces.
  # This does not matter mathematically, but it changes the order of floating
  # point summation in later correlations (~1e-17), which the search over
  # intermediate correlation matrices can amplify when eigenvalues are nearly
  # tied (e.g., symmetric blocks); keeping it identical keeps results identical
  build <- function(cor_matrix, standardize = FALSE){

    # Shared loadings (bounded, positive first loading)
    shared_load <- .cd_factor_analysis(cor_matrix, n_factors)
    shared_load[shared_load > 1] <- 1
    shared_load[shared_load < -1] <- -1
    if(shared_load[1, 1] < 0){
      shared_load <- shared_load * -1
    }

    # Unique loadings
    h2 <- rowSums(shared_load * shared_load)
    unique_load <- sqrt(ifelse(h2 < 1, 1 - h2, 0))

    # Normally distributed data
    normal <- tcrossprod(shared_comp, shared_load) +
      unique_comp * rep(unique_load, each = N)
    if(standardize){
      normal <- scale(normal)
    }

    # Impose observed distributions: for each variable in turn, order the
    # rows by that variable and assign its sorted distribution (the original
    # re-sorts the whole matrix each time; here only the row order is tracked)
    order_rows <- seq_len(N)
    imposed <- normal
    for(i in seq_len(k)){
      order_rows <- order_rows[sort.list(normal[order_rows, i])]
      imposed[order_rows, i] <- distributions[, i]
    }

    # Return population (rows in the original implementation's order)
    return(imposed[order_rows, , drop = FALSE])

  }

  # Search for the intermediate correlation matrix
  intermediate_cor <- target_cor
  while(trials_without_improvement < max_trials){

    # Population from current intermediate correlations
    population <- build(intermediate_cor)

    # Residual correlations
    reproduced_cor <- .cd_cor(population, cor_type)
    residual_cor <- target_cor - reproduced_cor
    rmsr <- sqrt(sum(residual_cor[lower]^2) / n_lower)

    # Update
    if(rmsr < best_rmsr){
      best_rmsr <- rmsr
      best_cor <- intermediate_cor
      best_res <- residual_cor
      intermediate_cor <- intermediate_cor + initial_multiplier * residual_cor
      trials_without_improvement <- 0
    }else{
      trials_without_improvement <- trials_without_improvement + 1
      intermediate_cor <- best_cor +
        initial_multiplier * 0.5^trials_without_improvement * best_res
    }

  }

  # Final population (standardized as in the original)
  return(build(best_cor, standardize = TRUE))

}

# Comparison data approach (Ruscio & Roche, 2012)
#
# Data = N x k data matrix; F.Max = largest number of factors considered;
# N.Pop = population size; N.Samples = samples per population;
# Alpha = significance level for adding a factor; Spearman = rank correlations;
# use = handling of missing values for the observed correlations
comparison_data <- function(
    Data, F.Max = 8, N.Pop = 10000, N.Samples = 500, Alpha = 0.30,
    Spearman = FALSE, use = "pairwise.complete.obs"
)
{

  # Dimensions
  Data <- as.matrix(Data)
  N <- nrow(Data)
  k <- ncol(Data)
  cor_type <- if(Spearman) "spearman" else "pearson"

  # Quantities that do not depend on the number of factors
  target_cor <- stats::cor(Data, method = cor_type, use = use)
  eigs_data <- eigen(target_cor, symmetric = TRUE, only.values = TRUE)$values
  distribution_values <- lapply(seq_len(k), function(i){na.omit(Data[, i])})

  # Search over number of factors
  rmsr_eigs <- matrix(0, nrow = N.Samples, ncol = F.Max)
  significant <- TRUE
  factors <- 1

  while((factors <= F.Max) && significant){

    # Comparison population
    population <- .cd_generate(
      Data, n_factors = factors, N = N.Pop, cor_type = cor_type,
      target_cor = target_cor, distribution_values = distribution_values
    )

    # Fit of eigenvalues from samples of the population
    for(j in seq_len(N.Samples)){
      sample_cor <- .cd_cor(
        population[sample.int(N.Pop, size = N, replace = TRUE), , drop = FALSE],
        cor_type
      )
      eigs_sample <- eigen(sample_cor, symmetric = TRUE, only.values = TRUE)$values
      rmsr_eigs[j, factors] <- sqrt(sum((eigs_sample - eigs_data)^2) / k)
    }

    # Is the added factor a significant improvement?
    if(factors > 1){
      significant <- stats::wilcox.test(
        rmsr_eigs[, factors], rmsr_eigs[, factors - 1], "less"
      )$p.value < Alpha
    }
    if(significant){
      factors <- factors + 1
    }

  }

  # Return number of factors
  return(factors - 1)

}

# Compile
comparison_data <- cmpfun(comparison_data)

#### Factor Forest ####
#
# Estimates the number of dimensions (factors) with the pre-trained xgboost
# model of Goretzko & Buhner (2020; https://osf.io/mvrau/), applied to ordinal
# data as in Goretzko & Buhner (2022; Applied Psychological Measurement).
#
# Feature extraction follows `analyze()` in `factor_forest/Evaluation.R` (the
# authors' own evaluation code for ordinal data). The model's 184 features are:
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
# The model is `factor_forest/tunedxgb.json`, a one-time conversion of
# `tunedxgb.rds` that current xgboost versions can read (see
# `factor_forest/convert_model.md`).
#
# Usage (inside `obtain_methods()`, re-using results that are already estimated):
#   factor_forest(simulated$data, R = R, pa = pa, kaiser = kaiser, ekc = ekc, cd = cd)
# Anything not supplied is computed internally, so `factor_forest(data)` also works.

# Folder this script was sourced from (NULL when it cannot be determined)
.methods_directory <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e){NULL})

# Cache so the model is only read once per R process
.factor_forest_cache <- new.env(parent = emptyenv())

# Load (and cache) the pre-trained model
#
# Searches (in order): `options(factor_forest.dir = "<folder>")`, the folder
# next to this script, and the project-root convention `./source/factor_forest`
factor_forest_model <- function(path = NULL)
{

  # Locate the model
  if(is.null(path)){

    # Candidate folders
    candidates <- file.path(
      c(
        getOption("factor_forest.dir"),
        if(!is.null(.methods_directory)){file.path(.methods_directory, "factor_forest")},
        "./source/factor_forest", "./factor_forest", "."
      ),
      "tunedxgb.json"
    )

    # Use the first that exists
    path <- Filter(file.exists, candidates)[1]
    if(is.na(path)){
      stop(
        "Could not find `tunedxgb.json`. Set `options(factor_forest.dir = \"<folder>\")`",
        call. = FALSE
      )
    }

  }

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

# Check that a supplied result is a single number
.factor_forest_check <- function(x, name)
{

  # Stop when not a single (non-missing) number
  if(!is.null(x) && !(is.numeric(x) && length(x) == 1 && !is.na(x))){
    stop(
      "`", name, "` must be a single number (the result of the corresponding method)",
      call. = FALSE
    )
  }

}

# Compute the model's features (one-row data frame)
#
# `data`: N x p matrix or data frame of item responses
# `R`: optional p x p correlation matrix (e.g., polychoric already computed);
#      computed from `data` when NULL
# `pa`: optional result of `parallel_analysis()`; re-used only when it includes
#       the "paf_fa" variant (the default `type` includes it), otherwise that
#       variant is computed
# `kaiser`, `ekc`, `cd`: optional results of `kaiser_guttman()`,
#       `empirical_kaiser()`, and `comparison_data()`; computed when NULL
#       (`cd` is by far the most expensive and is the most useful to supply)
# `cor`: "poly" (ordinal data; default) or "pearson" (continuous data)
factor_forest_features <- function(
    data, R = NULL, pa = NULL, kaiser = NULL, ekc = NULL, cd = NULL,
    cor = c("poly", "pearson"), use = "pairwise.complete.obs", F.Max = 8
)
{

  # Check supplied results
  .factor_forest_check(kaiser, "kaiser")
  .factor_forest_check(ekc, "ekc")
  .factor_forest_check(cd, "cd")

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
  if(is.null(ekc)){
    ekc <- empirical_kaiser(N = N, eigenvalues = eigval)
  }

  # Eigenvalue based features (including the Kaiser-Guttman rule)
  if(is.null(kaiser)){
    kaiser <- kaiser_guttman(eigenvalues = eigval)
  }
  eiggreater1 <- kaiser
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
  #
  # A supplied result is only usable if it includes the "paf_fa" variant
  if(is.null(pa) || !("paf_fa" %in% names(pa))){
    pa <- parallel_analysis(
      data, corr = switch(cor, "poly" = "polychoric", "pearson" = "pearson"),
      R = dat_cor, type = "paf_fa"
    )
  }
  pa_solution <- pa[["paf_fa"]]
  fa_eigval <- attr(pa, "eigenvalues")$paf_fa
  if(length(fa_eigval) != p){
    stop(
      "`pa` has ", length(fa_eigval), " factor analysis eigenvalues but the data ",
      "have ", p, " variables", call. = FALSE
    )
  }

  # Pad (or truncate) eigenvalues to the model's expected length
  eigval <- .factor_forest_pad(eigval)
  fa_eigval <- .factor_forest_pad(fa_eigval)
  names(eigval) <- paste0("eigval", seq_along(eigval))
  names(fa_eigval) <- paste0("fa_eigval", seq_along(fa_eigval))

  # Comparison data approach
  if(is.null(cd)){
    cd <- comparison_data(Data = data, F.Max = F.Max, use = use)
  }

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

# Compile
factor_forest_features <- cmpfun(factor_forest_features)

# Estimate the number of factors
#
# Returns an integer between 1 and 8 (the model's class labels)
#
# Arguments `R`, `pa`, `kaiser`, `ekc`, and `cd` are optional, already estimated
# results (see `factor_forest_features()`); `model` is a pre-trained booster
factor_forest <- function(
    data, R = NULL, pa = NULL, kaiser = NULL, ekc = NULL, cd = NULL,
    model = NULL, cor = c("poly", "pearson"),
    use = "pairwise.complete.obs", F.Max = 8
)
{

  # Obtain (cached) pre-trained model
  if(is.null(model)){
    model <- factor_forest_model()
  }

  # Compute features
  features <- factor_forest_features(
    data, R = R, pa = pa, kaiser = kaiser, ekc = ekc, cd = cd,
    cor = cor, use = use, F.Max = F.Max
  )

  # Predict class probabilities (uses the best early-stopping iteration)
  probabilities <- matrix(
    suppressWarnings(predict(model, as.matrix(features))),
    nrow = 1
  )

  # Return number of factors (class labels are 1, ..., 8)
  return(max.col(probabilities, ties.method = "first"))

}

# Compile
factor_forest <- cmpfun(factor_forest)

# Obtain methods
obtain_methods <- function(simulated)
{

  # Get variable names
  variables <- colnames(simulated$data)

  # Get sample size
  N <- dim(simulated$data)[1]

  # Estimate correlations
  R <- polychoric(simulated$data)

  # EGA
  ega <- EGA(R, n = N, plot.EGA = FALSE)

  # Initialize unidimensional flag
  unidimensional <- ega$n.dim == 1

  # Check for flag
  uni_solution <- rep(1, length(variables))

  ## Walktrap
  walktrap <- EGAnet:::swiftelse(
    unidimensional, uni_solution,
    community.detection(ega$network, algorithm = "walktrap")
  )

  ## Louvain
  louvain <- EGAnet:::swiftelse(
    unidimensional, uni_solution,
    community.consensus(ega$network)
  )

  ## PRISM
  prism <- EGAnet:::swiftelse(
    unidimensional, uni_solution,
    community.prism(ega$network)
  )

  # Estimate parallel analysis
  pa <- parallel_analysis(simulated$data, corr = "polychoric")

  # Estimate Kaiser-Guttman
  kaiser <- kaiser_guttman(R = R)

  # Estimate empirical Kaiser
  ekc <- empirical_kaiser(R = R, N = N)

  # Estimate comparison data
  cd <- comparison_data(Data = simulated$data)

  # Estimate factor forest
  ff <- factor_forest(
    simulated$data, R = R, pa = pa, kaiser = kaiser, ekc = ekc, cd = cd
  )

  # Obtain results
  results <- list(
    walktrap = walktrap, louvain = louvain, prism = prism,
    pca = pa[["pca"]], paf = pa[["paf"]], pafm = pa[["paf_fa"]],
    kaiser = kaiser, ekc = ekc, cd = cd, ff = ff
  )

  # Add correlations
  attributes(results)$R <- R

  # Return results
  return(results)

}

# Compile
obtain_methods <- cmpfun(obtain_methods)
