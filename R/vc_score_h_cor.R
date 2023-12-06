#'Computes variance component score test statistic for homogeneous effects
#'
#'This function computes the variance component score test statistics assuming homogeneous
#'effects of the testing variable across all genes in the set
#'
#'@keywords internal
#'
#'@param y a numeric matrix of dim \code{g x n} containing the raw or
#'normalized RNA-seq counts for g genes from \code{n} samples.
#'
#'@param x a numeric design matrix of dim \code{n x p} containing the \code{p}
#'covariates to be adjusted for. First column should contain 1s if an intercept is desired.
#'
#'@param indiv a vector of length \code{n} containing the information for
#'attributing each sample to one of the studied individuals. Coerced
#'to be a \code{factor}.
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the
#'\code{K} longitudinal variables to be tested (typically a vector of time
#'points or functions of time).
#'
#'#'@param use_phi a logical flag indicating whether conditional means should be
#'conditioned on \code{phi} and on covariate(s) \code{x}, or on \code{x} alone.
#'Default is \code{TRUE} in which case conditional means are estimated
#'conditionally on both \code{x} and \code{phi}.
#'
#'#'@param gene_based a logical flag indicating if weights have been estimated
#'at the gene-level. Default is \code{FALSE}, indicating the input weights 
#'are observation level. In the case where weights are gene-specific and not
#'observation specific, setting this to \code{TRUE} will result in more efficient 
#'computation.
#'
#'@param Sigma a list of \code{n} matrices of size \code{g x g} giving the within-set
#'gene-gene correlation structure and heteroskedasticity weights (computed by
#'the function GS_cor()).
#'
#'@param Sigma_xi a matrix of size \code{K x K} containing the covariance matrix
#'of the \code{K} random effects corresponding to \code{phi}.
#'
#'@param na_rm logical: should missing values (including \code{NA} and
#'\code{NaN}) be omitted from the calculations? Default is \code{FALSE}.
#'
#'@return A list with the following elements:\itemize{
#'   \item \code{score}: approximation of the set observed score
#'   \item \code{q}: observation-level contributions to the score
#'   \item \code{gene_scores}: approximation of the individual gene scores
#' }
#'
#'
#'@export
vc_score_h_cor <- function(y, x, indiv, phi, use_phi = T, Sigma, gene_based = F, Sigma_xi = diag(ncol(phi)),
                       na_rm = FALSE) {

    ## validity checks
    if ((sum(!is.finite(unlist(Sigma)))) > 0) {
        stop("At least 1 non-finite covariance entry in Sigma")
    }
  
  # checking weights are truly as user specified:
  if(gene_based == TRUE & !all(Sigma[[2]] == Sigma[[1]])){
    warning("You indicated weights to be calculated on the gene level, however they appear",
            " to be on the observation level. Setting gene_based = FALSE.")
    gene_based = FALSE
  }

  stopifnot(is.matrix(y))
  stopifnot(is.matrix(x))
  stopifnot(is.matrix(phi))
  
  g <- nrow(y)  # the number of genes measured
  n <- ncol(y)  # the number of samples measured
  qq <- ncol(x)  # the number of covariates
  n_t <- ncol(phi)  # the number of test variables
  stopifnot(nrow(x) == n)
  stopifnot(length(Sigma) == n)
  stopifnot(dim(Sigma[[1]]) == g)
  stopifnot(nrow(phi) == n)
  stopifnot(length(indiv) == n)
  
  # the number of random effects
  if (length(Sigma_xi) == 1) {
    K <- 1
    Sigma_xi <- matrix(Sigma_xi, K, K)
  } else {
    K <- nrow(Sigma_xi)
    stopifnot(ncol(Sigma_xi) == K)
  }
  stopifnot(n_t == K)
  
  ## data formating ------
  indiv <- as.factor(indiv)
  nb_indiv <- length(levels(indiv))
  
  y_T <- t(y)
  
  # Fitting generalised OLS model to individual samples of the g genes in the set
  # with covariance of residual errors given by Sigma
  # effect of covariates (and possibly phi)
  xphi <- if (use_phi) {
    cbind(x, phi)
  } else {
    x
  }
  
  y_mu = matrix(ncol = n, nrow = g) # Matrix to hold centralised expression
  q_individual = matrix(ncol = g, nrow = n) # Matrix to hold observation contributions to stat
  for (i in c(1:n)){
    X_diag = t(as.matrix(bdiag(rep(list(xphi[i,]), g)))) #Block diagonal design matrix
    y_i = y[,i]
    sigma_i = Sigma[[i]] # ith individual's covariance matrix
    sigma_inv = Ginv(sigma_i) # (generalised) inverse of covariance matrix
    xt_sig_x_inv = Ginv((t(X_diag) %*% sigma_i %*% X_diag))
    y_mu[,i] = y_i - X_diag %*% xt_sig_x_inv %*% t(X_diag) %*% sigma_inv %*% y_i 
    #gene expression minus generalised OLS estimate
    phi_diag = t(as.matrix(bdiag(rep(list(phi[i]), g)))) # Phi block diagonal design matrix
    q_individual[i,] = t(y_mu[,i]) %*% sigma_inv %*% phi_diag # observation level contributions
  }
  
  q <- (n*g)^{-1/2} *colSums(q_individual, na.rm = na_rm) #gene level scores
  Q = t(q) %*% q   # set score
  
  
  return(list(score = Q, q = q_individual, gene_scores = q))
}
