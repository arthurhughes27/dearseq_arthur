#'Estimation of gene-set correlation structure
#'
#'Computes estimated intra-gene set inter-gene correlation structure incorporating 
#'heteroskedasticity weights 
#'
#'@param y a numeric matrix of size \code{p x n} containing the raw RNA-seq
#'counts or preprocessed expression from \code{n} samples for \code{p} genes in 
#'a singular gene set .
#'
#'@param x a numeric matrix of size \code{n x q} containing the model
#'covariate(s) from \code{n} samples (design matrix).
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the K
#'variable(s) whose effects on the expression inside the gene set 
#'are of interest( e.g. bases of time).
#'
#'@param indiv a vector of length \code{n} containing the information for
#'attributing each sample to one of the studied individuals. Coerced
#'to be a \code{factor}. Default is \code{NULL}, in which case it is assumed that
#'each row represents a distinct individual.
#'
#'@param w: a matrix \code{n x p} containing the computed precision 
#'weights. Default is \code{NULL}, in which case heteroskedasticity is 
#'not taken into account. 
#'
#'@param use_phi a logical flag indicating whether conditional means should be
#'conditioned on \code{phi} and on covariate(s) \code{x}, or on \code{x} alone.
#'Default is \code{TRUE} in which case conditional means are estimated
#'conditionally on both \code{x} and \code{phi}.
#'
#'@param preprocessed a logical flag indicating whether the expression data have
#'already been preprocessed (e.g. log2 transformed). Default is \code{FALSE}, in
#'which case \code{y} is assumed to contain raw counts and is normalized into
#'log(counts) per million.
#'
#'@param verbose a logical flag indicating whether informative messages are
#'printed during the computation. Default is \code{TRUE}.
#'
#'@param na.rm logical: should missing values (including \code{NA} and
#'\code{NaN}) be omitted from the calculations? Default is \code{FALSE}.
#'
#'@return \code{Sigma_list}: a list of \code{n} matrices of size \code{p x p} 
#' containing the computed observation level covariance matrices of the residual 
#' errors conditional on the covariates (and possibly \code{phi}).
#'
#'@author Arthur Hughes
#'
#'@seealso \code{\link[stats]{bandwidth}} \code{\link{density}}
#'
#'@import ggplot2
#'@import matlib
#'@importFrom stats bw.bcv bw.nrd0 bw.nrd bw.SJ bw.ucv dnorm approx sd pnorm
#'qnorm na.omit
#'@importFrom matlib Ginv
#'@export

GS_cor <- function(y, x, phi = NULL, indiv, w = NULL, use_phi = TRUE, preprocessed = FALSE, 
                   verbose = TRUE, na.rm = FALSE){
  
  ## dimensions & validity checks
  stopifnot(is.matrix(y)) #y, x, and phi must be matrices
  stopifnot(is.matrix(x))
  stopifnot(is.null(phi) | is.matrix(phi)) 
  
  p <- nrow(y)  # the number of genes in the set
  n <- ncol(y)  # the number of samples measured
  q <- ncol(x)  # the number of covariates
  K <- ncol(phi) # the number of test variables 
  
  stopifnot(nrow(x) == n) # x must have samples as rows
  stopifnot(length(indiv) == n) # vector of sample indices must have n elements
  if(use_phi){
    stopifnot(nrow(phi) == n) # phi must have samples as rows 
  }
  
  if (is.null(w)){
    w = matrix(1, nrow = p, ncol = n) # if no weights specified, heteroskedasticity not accounted for
  } else {
    stopifnot(nrow(w) == p | ncol(w) == n) #else weights must have column samples and row genes
  }
  
  indiv <- as.factor(as.numeric(as.factor(indiv)))  #coercing sample index vector to factors
  n_indiv <- length(levels(indiv)) #number of individuals under study
  n_i = (summary(indiv, maxsum = n_indiv) %>% as.vector()) #number of observations per individual
  t = max(n_i) # max number of observations 

  # removing genes never observed:
  observed <- which(rowSums(y, na.rm = TRUE) != 0) 
  nb_g_sum0 <- length(observed) - p
  if (nb_g_sum0 > 0) {
    warning(nb_g_sum0, " y rows sum to 0 (i.e. are never observed)",
            "and have been removed")
  }
  
  #if counts are not already normalised, convert to log-counts per million
  if (preprocessed) {
    y_lcpm <- t(y[observed, ])
  } else {
    y_lcpm <- t(apply(y[observed, ], MARGIN = 2, function(v) {
      log2((v + 0.5)/(sum(v) + 1) * 10^6)
    }))
  }
  p <- ncol(y_lcpm) # number of genes in set post null removal 
  
  # Fitting genewise OLS model to estimate gene expression without mean 
  # effect of covariates (and possibly phi)
  xphi <- if (use_phi) {
    cbind(x, phi)
  } else {
    x
  }
  
  if (na.rm) {
    y_lcpm0 <- y_lcpm
    y_lcpm0[is.na(y_lcpm0)] <- 0 # remove na counts if option selected
    alpha_ols <- solve(crossprod(xphi)) %*% t(xphi) %*% y_lcpm0 # least squares estimate of beta
    R_gene = y_lcpm0 - xphi %*% alpha_ols # residual expression after mean effect removed
  } else {
    alpha_ols <- solve(crossprod(xphi)) %*% t(xphi) %*% y_lcpm # least squares estimate of beta
    R_gene = y_lcpm - xphi %*% alpha_ols # residuals w.r.t. genes
  }
  
  #residuals w.r.t. obs
  R_obs = matrix(NA,nrow = n_indiv*p, ncol = t) # initialise
  for (g in c(1:p)){
    for (i in c(1:n_indiv)){
      R_obs[(g-1)*n_indiv+i,c(1:n_i[i])] = R_gene[indiv == i,g]
    } 
  }
  
  cor_gene = cov(R_gene) # gene-wise correlation structure
  cor_obs = cor(R_obs, use = "pairwise.complete.obs") # observation-wise correlation structure
  
  # Plug in estimator of correlation matrix
  P = matrix(0, nrow = t*p, ncol = t*p) # empty pT*pT matrix
  for (g1 in c(1:p)){
    for (g2 in c(1:p)){
      P[c(((g1-1)*t+1):(g1*t)),
        c(((g2-1)*t+1):(g2*t))] = cor_gene[g1,g2]*cor_obs 
      # each block of P is the observation-wise matrix multiplied 
      # by a single element of the genewise matrix
    }
  }

  #Plug in estimator of covariance matrix incorporating weights
  sigma_list = vector("list", n_indiv) #store individual covariance matrices
  for (i in c(1:n_indiv)){ #for all individuals
    wcol = matrix(0,nrow = p*n_i[i])
    for (g in c(1:p)){
      wcol[c(((g-1)*n_i[i]+1):(g*n_i[i]))] = w[g,indiv == i]
    }
    
    sigma_list[[i]] = matrix(0, ncol = p*n_i[i], nrow = p*n_i[i]) # initialise matrices
    for (j in c(1:(p*n_i[i]))){
      sigma_list[[i]][,j] = P[c(1:(n_i[i]*p)),j]*wcol[j]^{-1/2}*wcol^{-1/2} # ith covariance matrix
    }
  }
  return(sigma_list)
}


