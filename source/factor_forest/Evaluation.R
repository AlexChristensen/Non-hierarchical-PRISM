library(batchtools)
library(psych)
library(mvtnorm)
library(ineq)
library(purrr)
library(mlr)
library(data.table)
library(BBmisc)

reg <- makeExperimentRegistry("evaluation_factor_forest", seed = 95558)

# load packages on each node ("slaves")

reg$packages <- c("psych","mvtnorm","ineq","purrr","mlr",
            "data.table", "BBmisc")

# source xgboost functions for specific "early-stopping" implementation

source("xgb-functions.R")

# function to create loading matrix 

load_mat <- function(k, p, vpf, pf, sf){

  switch(pf,
         "low" = {
           pf_low = 0.35 
           pf_upper = 0.5
           },
         "medium" = {
           pf_low = 0.5
           pf_upper = 0.65
         },
         "high" = {
           pf_low = 0.65
           pf_upper = 0.8
         })
  switch(sf,
         "none" = {
           sf_low = 0 
           sf_upper = 0
         },
         "low" = {
           sf_low = 0
           sf_upper = 0.1
         },
         "medium" = {
           sf_low = 0.1
           sf_upper = 0.2
         })
  
  x <- runif(p, pf_low, pf_upper) 
  y <- runif(p*(k-1) , sf_low, sf_upper)
  
  i <- 1:(p)

  j <- rep(1:k, each=vpf) 
  
  L <- matrix(NA, p, k) 
  L[cbind(i, j)] <- x 
  L[is.na(L)] <- y
  L
}


# CD approach after Ruscio and Roche 2012 (internal seed deleted)

EFA.Comp.Data <- function(Data, F.Max, N.Pop = 10000, N.Samples = 500, Alpha = .30, Graph = F, Spearman = F)
{
  # Data = N (sample size) x k (number of variables) data matrix
  # F.Max = largest number of factors to consider
  # N.Pop = size of finite populations of comparison data (default = 10,000 cases)
  # N.Samples = number of samples drawn from each population (default = 500)
  # Alpha = alpha level when testing statistical significance of improvement with add'l factor (default = .30) 
  # Graph = whether to plot the fit of eigenvalues to those for comparison data (default = F)
  # Spearman = whether to use Spearman rank-order correlations rather than Pearson correlations (default = F)
  
  N <- dim(Data)[1]
  k <- dim(Data)[2]
  if (Spearman) Cor.Type <- "spearman" else Cor.Type <- "pearson"
  cor.Data <- cor(Data, method = Cor.Type)
  Eigs.Data <- eigen(cor.Data)$values
  RMSR.Eigs <- matrix(0, nrow = N.Samples, ncol = F.Max)
  Sig <- T
  F.CD <- 1
  while ((F.CD <= F.Max) & (Sig))
  {
    Pop <- GenData(Data, N.Factors = F.CD, N = N.Pop, Cor.Type = Cor.Type)
    for (j in 1:N.Samples)
    {
      Samp <- Pop[sample(1:N.Pop, size = N, replace = T),]
      cor.Samp <- cor(Samp, method = Cor.Type)
      Eigs.Samp <- eigen(cor.Samp)$values
      RMSR.Eigs[j,F.CD] <- sqrt(sum((Eigs.Samp - Eigs.Data) * (Eigs.Samp - Eigs.Data)) / k)
    }
    if (F.CD > 1) Sig <- (wilcox.test(RMSR.Eigs[,F.CD], RMSR.Eigs[,(F.CD - 1)], "less")$p.value < Alpha)
    if (Sig) F.CD <- F.CD + 1
  }
  if (Graph)
  {
    if (Sig) x.max <- F.CD - 1
    else x.max <- F.CD
    ys <- apply(RMSR.Eigs[,1:x.max], 2, mean)
    plot(x = 1:x.max, y = ys, ylim = c(0, max(ys)), xlab = "Factor", ylab = "RMSR Eigenvalue", type = "b", 
         main = "Fit to Comparison Data")
    abline(v = F.CD - 1, lty = 3)
  }
  return(F.CD - 1)
}


################################################################################################################
GenData <- function(Supplied.Data, N.Factors, N, Max.Trials = 5, Initial.Multiplier = 1, Cor.Type){
  
  k <- dim(Supplied.Data)[2]
  Data <- matrix(0, nrow = N, ncol = k)            # Matrix to store the simulated data
  Distributions <- matrix(0, nrow = N, ncol = k)   # Matrix to store each variable's score distribution
  Iteration <- 0                                   # Iteration counter
  Best.RMSR <- 1                                   # Lowest RMSR correlation
  Trials.Without.Improvement <- 0                  # Trial counter
  
  # Generate distribution for each variable (step 2) -------------------------------------------------------------
  
  for (i in 1:k)
    Distributions[,i] <- sort(sample(Supplied.Data[,i], size = N, replace = T))
  
  # Calculate and store a copy of the target correlation matrix (step 3) -----------------------------------------
  
  Target.Corr <- cor(Supplied.Data, method = Cor.Type)
  Intermediate.Corr <- Target.Corr
  
  # Generate random normal data for shared and unique components, initialize factor loadings (steps 5, 6) --------
  
  Shared.Comp <- matrix(rnorm(N * N.Factors, 0, 1), nrow = N, ncol = N.Factors)
  Unique.Comp <- matrix(rnorm(N * k, 0, 1), nrow = N, ncol = k)
  Shared.Load <- matrix(0, nrow = k, ncol = N.Factors)
  Unique.Load <- matrix(0, nrow = k, ncol = 1)
  
  # Begin loop that ends when specified number of iterations pass without improvement in RMSR correlation --------
  
  while (Trials.Without.Improvement < Max.Trials)
  {
    Iteration <- Iteration + 1
    
    # Calculate factor loadings and apply to reproduce desired correlations (steps 7, 8) ---------------------------
    
    Fact.Anal <- Factor.Analysis(Intermediate.Corr, Corr.Matrix = TRUE, N.Factors = N.Factors, Cor.Type = Cor.Type)
    if (N.Factors == 1) Shared.Load[,1] <- Fact.Anal$loadings
    else 
      for (i in 1:N.Factors)
        Shared.Load[,i] <- Fact.Anal$loadings[,i]
      Shared.Load[Shared.Load > 1] <- 1
      Shared.Load[Shared.Load < -1] <- -1
      if (Shared.Load[1,1] < 0) Shared.Load <- Shared.Load * -1
      for (i in 1:k)
        if (sum(Shared.Load[i,] * Shared.Load[i,]) < 1) Unique.Load[i,1] <- 
        (1 - sum(Shared.Load[i,] * Shared.Load[i,]))
      else Unique.Load[i,1] <- 0
      Unique.Load <- sqrt(Unique.Load)
      for (i in 1:k)
        Data[,i] <- (Shared.Comp %*% t(Shared.Load))[,i] + Unique.Comp[,i] * Unique.Load[i,1]
      
      # Replace normal with nonnormal distributions (step 9) ---------------------------------------------------------
      
      for (i in 1:k)
      {
        Data <- Data[sort.list(Data[,i]),]
        Data[,i] <- Distributions[,i]
      }
      
      # Calculate RMSR correlation, compare to lowest value, take appropriate action (steps 10, 11, 12) --------------
      
      Reproduced.Corr <- cor(Data, method = Cor.Type)
      Residual.Corr <- Target.Corr - Reproduced.Corr
      RMSR <- sqrt(sum(Residual.Corr[lower.tri(Residual.Corr)] * Residual.Corr[lower.tri(Residual.Corr)]) / 
                     (.5 * (k * k - k)))
      if (RMSR < Best.RMSR)
      {
        Best.RMSR <- RMSR
        Best.Corr <- Intermediate.Corr
        Best.Res <- Residual.Corr
        Intermediate.Corr <- Intermediate.Corr + Initial.Multiplier * Residual.Corr
        Trials.Without.Improvement <- 0
      }
      else 
      {
        Trials.Without.Improvement <- Trials.Without.Improvement + 1
        Current.Multiplier <- Initial.Multiplier * .5 ^ Trials.Without.Improvement
        Intermediate.Corr <- Best.Corr + Current.Multiplier * Best.Res
      }
  }
  
  Fact.Anal <- Factor.Analysis(Best.Corr, Corr.Matrix = TRUE, N.Factors = N.Factors, Cor.Type = Cor.Type)
  if (N.Factors == 1) Shared.Load[,1] <- Fact.Anal$loadings
  else
    for (i in 1:N.Factors)
      Shared.Load[,i] <- Fact.Anal$loadings[,i]
  Shared.Load[Shared.Load > 1] <- 1
  Shared.Load[Shared.Load < -1] <- -1
  if (Shared.Load[1,1] < 0) Shared.Load <- Shared.Load * -1
  for (i in 1:k)
    if (sum(Shared.Load[i,] * Shared.Load[i,]) < 1) Unique.Load[i,1] <-
    (1 - sum(Shared.Load[i,] * Shared.Load[i,]))
  else Unique.Load[i,1] <- 0
  Unique.Load <- sqrt(Unique.Load)
  for (i in 1:k)
    Data[,i] <- (Shared.Comp %*% t(Shared.Load))[,i] + Unique.Comp[,i] * Unique.Load[i,1]
  Data <- apply(Data, 2, scale) # standardizes each variable in the matrix
  for (i in 1:k)
  {
    Data <- Data[sort.list(Data[,i]),]
    Data[,i] <- Distributions[,i]
  }
  
  return(Data)
}

################################################################################################################
Factor.Analysis <- function(Data, Corr.Matrix = FALSE, Max.Iter = 50, N.Factors = 0, Cor.Type)
{
  Data <- as.matrix(Data)
  k <- dim(Data)[2]
  if (N.Factors == 0)
  {
    N.Factors <- k
    Determine <- T
  }
  else Determine <- F
  if (!Corr.Matrix) Cor.Matrix <- cor(Data, method = Cor.Type)
  else Cor.Matrix <- Data
  Criterion <- .001
  Old.H2 <- rep(99, k)
  H2 <- rep(0, k)
  Change <- 1
  Iter <- 0
  Factor.Loadings <- matrix(nrow = k, ncol = N.Factors)
  while ((Change >= Criterion) & (Iter < Max.Iter))
  {
    Iter <- Iter + 1
    Eig <- eigen(Cor.Matrix)
    L <- sqrt(Eig$values[1:N.Factors])
    for (i in 1:N.Factors)
      Factor.Loadings[,i] <- Eig$vectors[,i] * L[i]
    for (i in 1:k)
      H2[i] <- sum(Factor.Loadings[i,] * Factor.Loadings[i,])
    Change <- max(abs(Old.H2 - H2))
    Old.H2 <- H2
    diag(Cor.Matrix) <- H2
  }
  if (Determine) N.Factors <- sum(Eig$values > 1)
  return(list(loadings = Factor.Loadings[,1:N.Factors], factors = N.Factors))
}


# source pre-trained model

modxgb <- readRDS(file = "tunedxgb.rds")

## export support functions (functions have to be available on each node)

batchExport(export = list(load_mat = load_mat, EFA.Comp.Data = EFA.Comp.Data,
                          GenData = GenData, Factor.Analysis = Factor.Analysis,
                          modxgb = modxgb,
                          makeRLearner.classif.xgboost.earlystop = makeRLearner.classif.xgboost.earlystop,
                          trainLearner.classif.xgboost.earlystop = trainLearner.classif.xgboost.earlystop,
                          predictLearner.classif.xgboost.earlystop = predictLearner.classif.xgboost.earlystop,
                          createDMatrixFromTask = createDMatrixFromTask))

# simulate data

simulation <- function(data, job,  # mandatory arguments
                       N,vpf, k, rho, pf, sf, distr, cat){
  # simulate data

  # simulate loading pattern
  
  p = vpf*k
  L <- load_mat(k, p, vpf, pf, sf)
  fcor <- matrix(rho, k,k) + diag(1-rho,k,k)
  Sigma <- L%*%fcor%*%t(L) + diag(diag(diag(p)-L%*%fcor%*%t(L)))
  
  # simulate data with specific distribution depending on simulation condition
  
  dat <- switch(distr,
                "ordinal" = {
                    dat <- data.frame(mvtnorm::rmvnorm(N, mean = rep(0,p), sigma = Sigma))
                    dat <- data.frame(
                      lapply(dat,
                           FUN = function(x){cut(x, c(min(x)-0.01,
                                                      qnorm((1:(cat-1))/cat),
                                                      max(x) + 0.01), labels = FALSE,
                                                 right =FALSE)}))
                  dat
                },
                "ordinal_skewed" = {
                  dat <- data.frame(mvtnorm::rmvnorm(N, mean = rep(0,p), sigma = Sigma))
                  dat <- data.frame(
                    lapply(dat,
                           FUN = function(x){cut(x, c(min(x)-0.01,
                                                      sort(sample(c(-1,1),1)*qnorm((1:(cat-1))/(1.5*cat))),
                                                      max(x) + 0.01), labels = FALSE,
                                                 right =FALSE)}))
                  dat
                }
  )
 
 
  list(dat = dat, L = L, p = p, k = k, N = N, distr = distr)
}

addProblem("data_simulation", fun = simulation)

# add analyze function (major parts of this function are based on the script by
# Goretzko & Bühner, 2020 who provided the factor.forest() function in an OSF
# repository - https://osf.io/mvrau/)

analyze <- function(data, job, instance, ...){
  
  # use polychoric correlations instead of Pearson if data are ordinal
  
  dat_cor <- polychoric(instance$dat)$rho
  
  eigval <- eigen(dat_cor)$values
  vareig <- cumsum(eigval)/instance$p
  
  # empirical kaiser criterion
  
  lref <- rep(0,instance$p)
  for (i in 1:instance$p) {
    lref[i] <- max(((1 + sqrt(instance$p/instance$N))^2) * (instance$p-sum(lref))/(instance$p-i+1),1)
  }
  ekc <- which(eigval<=lref)[1]-1
  
  # calculate eigenvalue based features
  
  eiggreater1 <- sum(eigval > 1)  
  releig1 <- eigval[1]/instance$p
  releig2 <- sum(eigval[1:2])/instance$p
  releig3 <- sum(eigval[1:3])/instance$p
  eiggreater07 <- sum(eigval > 0.7)
  sdeigval <- sd(eigval)
  var50 <- min(which(vareig > 0.50))
  var75 <- min(which(vareig > 0.75))
  
  # calculate matrix norm features
  
  onenorm <- norm(dat_cor,"O")
  frobnorm <- norm(dat_cor,"F")
  maxnorm <- norm(dat_cor-diag(instance$p),"M")
  avgcor <- sum(abs(dat_cor-diag(instance$p)))/(instance$p*(instance$p-1))
  specnorm <- sqrt(eigen(t(dat_cor)%*%dat_cor)$values[1])
  
  smlcor <- sum(dat_cor <= 0.1)
  avgcom <- mean(smc(dat_cor))
  det <- det(dat_cor)
  
  # "inequality" features 
  
  KMO <- KMO(dat_cor)$MSA
  Gini <- ineq(lower.tri(dat_cor), type = "Gini")
  Kolm <- ineq(lower.tri(dat_cor), type = "Kolm")
  
  # parallel analysis
  
  pa <- fa.parallel(instance$dat,fa="fa",plot = FALSE, cor = "poly")
  pa_solution <- pa$nfact
  fa_eigval <- pa$fa.values
  
  # adding -1000 to shorter eigenvalue vectors 
  
  eigval[(length(eigval)+1):80] <- -1000
  fa_eigval[(length(fa_eigval)+1):80] <- -1000
  names(eigval) <- paste("eigval", 1:80, sep = "")
  names(fa_eigval) <- paste("fa_eigval", 1:80, sep = "")
  
  # comparison data approach (with maximum number of factors = 8)
  
  cd <- EFA.Comp.Data(Data = instance$dat, F.Max = 8)
  
  # combine features to data frame (input for trained model) 
  
  features <- cbind(data.frame(k = instance$k,N=instance$N,p=instance$p,eiggreater1,
                               releig1,releig2,releig3,eiggreater07,sdeigval,var50,
                               var75,onenorm,frobnorm,maxnorm, avgcor,specnorm, smlcor,
                               avgcom,det, KMO, Gini, Kolm, pa_solution, ekc, cd),
                               t(eigval), t(fa_eigval))
  
  # predict number of factors based on observed features
  
  pred_xgb <- predict(modxgb, newdata=features)
  
  # store results of parallel analysis, Kaiser-Guttman rule, factor forest,
  # empirical Kaiser criterion and comparison data approach)
  
  data.frame(pa_solution, guttman = eiggreater1,
             xgb = pred_xgb$data$response, ekc, cd)
}

addAlgorithm("analyze", fun = analyze)

# create experimental conditions

cond <- CJ(N = c(250, 1000), vpf = c(4, 7), k = c(1,4,6),
           rho = c(0,0.3), pf = c("low", "medium","high"),
           sf = c("none","low","medium"),  distr = c("ordinal", "ordinal_skewed"),
           cat = 2:6)

simulation_design <- list(data_simulation = 
                            cond[!((cond$pf == "high" & cond$sf == "medium" & cond$k >= 6 & cond$rho >= 0.3) |
                                     (cond$k == 1 & cond$rho == 0.3) | (cond$k == 1 & sf != "none")),])

addExperiments(prob.designs = simulation_design, repls = 500)
