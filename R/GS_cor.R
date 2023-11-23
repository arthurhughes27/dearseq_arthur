#'Estimation of gene-set correlation structure
#'
#'Computes estimated intra-gene set inter-gene correlation structure incorporating 
#'heteroskedasticity weights 
#'
#'@param y a numeric matrix of size \code{G x n} containing the raw RNA-seq
#'counts or preprocessed expression from \code{n} samples for \code{G} genes in 
#'a singular gene set .
#'
#'@param x a numeric matrix of size \code{n x p} containing the model
#'covariate(s) from \code{n} samples (design matrix).
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the K
#'variable(s) whose effects on the expression inside the gene set 
#'are of interest( e.g. bases of time).
#'
#'@param weights: a matrix \code{n x G} containing the computed precision 
#'weights
#'
#'@param gene_based a logical flag indicating if weights have been estimated
#'at the gene-level. Default is \code{FALSE}, indicating the input weights 
#'are observation level. In the case where weights are gene-specific and not
#'observation specific, setting this to \code{TRUE} will result in more efficient 
#'computation.
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
#' errors conditional on the covariates (and possibly \code{phi}). In the case 
#' where gene_based = TRUE, returns a list of one matrix. 
#'
#'@author Arthur Hughes
#'
#'@seealso \code{\link[stats]{bandwidth}} \code{\link{density}}
#'
#'@import ggplot2
#'@importFrom stats bw.bcv bw.nrd0 bw.nrd bw.SJ bw.ucv dnorm approx sd pnorm
#'qnorm na.omit
#'@export

GS_cor <- function(y, x, phi = NULL, weights, gene_based = F, use_phi = TRUE, preprocessed = FALSE, verbose = TRUE,
                       na.rm = FALSE){
  
  ## dimensions & validity checks
  stopifnot(is.matrix(y)) #y, x, and phi must be matrices
  stopifnot(is.matrix(x))
  stopifnot(is.null(phi) | is.matrix(phi))
  g <- nrow(y)  # the number of genes in the set
  n <- ncol(y)  # the number of samples measured
  qq <- ncol(x)  # the number of covariates
  stopifnot(nrow(x) == n) 
  if(use_phi){
    stopifnot(nrow(phi) == n) # phi and x must have n measurements
  }
  
  # checking weights are truly as user specified:
  if(gene_based == TRUE & !all(apply(weights, 2, identical, weights[,1]))){
    warning("You indicated weights to be calculated on the gene level, however they appear",
            " to be on the observation level. Setting gene_based = FALSE.")
    gene_based = F
    }
  
  # removing genes never observed:
  observed <- which(rowSums(y, na.rm = TRUE) != 0)
  nb_g_sum0 <- length(observed) - g
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
  N <- length(y_lcpm)
  p <- ncol(y_lcpm)
  
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
    B_ols <- solve(crossprod(xphi)) %*% t(xphi) %*% y_lcpm0 # least squares estimate of beta
    R = y_lcpm0 - xphi %*% B_ols # residual expression after mean effect removed
  } else {
    B_ols <- solve(crossprod(xphi)) %*% t(xphi) %*% y_lcpm # least squares estimate of beta
    R = y_lcpm - xphi %*% B_ols # residual expression after mean effect removed
  }
  
  cor_R = cor(R) # correlation of residuals estimating the correlation structure in the set
  
  
  if (gene_based == F){
    sigma_list = vector("list", n) # list to store covariance matrices
    for (i in c(1:n)){ #for all individuals
      sigma_list[[i]] = matrix(nrow = p, ncol= p) #list of observation level covariance matrices
      for (j in c(1:p)){
        for (k in c(1:p)){
          sigma_list[[i]][j,k] = cor_R[j,k]*sqrt(weights[j,i])*sqrt(weights[k,i])
        # multiply sample correlation by square root of
        # individual level heteroskedasticity weights
        }
      }
    }
  } else {
    sigma_list = vector("list", 1) # if weights were gene level, only need to compute 1 matrix
    sigma_list[[1]] = matrix(nrow = p, ncol= p)
    for (j in c(1:p)){
      for (k in c(1:p)){
        sigma_list[[1]][j,k] = cor_R[j,k]*sqrt(weights[j,1])*sqrt(weights[k,1])
      }
    }
  }
  return(sigma_list)
}


