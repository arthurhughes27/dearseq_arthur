#'Computes variance component score test statistic for homogeneous trajectories
#'
#'This function computes the variance component score test statistics for a set taking into account
#'heteroskedasticity, intra-gene set correlation structures and intra-individual correlation 
#'in the case where we assume heterogeneous trajectories with respect to the covariates and testing
#'variables.
#'
#'@keywords internal
#'
#'@param y a numeric matrix of dim \code{p x n} containing the expression
#'for \code{p} genes from \code{n} samples.
#'
#'@param x a numeric design matrix of dim \code{n x q} containing the \code{q}
#'covariates to be adjusted for.
#'
#'@param indiv a vector of length \code{n} containing the information for
#'attributing each sample to one of the studied individuals. Coerced
#'to be a \code{factor}. By default, it is assumed that each row represents a distinct individual.
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the
#'\code{K} variables to be tested.
#'
#'@param w a vector of length \code{n} containing the weights for the \code{n}
#'samples. Default is \code{NULL}, in which case heteroskedasticity is 
#'not taken into account.
#'
#'@param Sigma a list of matrices giving the plug-in estimate of the intra-set and 
#'inter-observational covariance structures, incorporating the heteroskedasticity weights
#'(computed by the function GS_cor()). The list should be the same length as the number of 
#'distinct individuals under study. Each element of the list corresponds to individual \code{i}'s 
#'covariance matrix, and has dimension \code{n_{i} x n_{i}}, where \code{n_{i}} is the number 
#'of samples corresponding to individual \code{i}. Default is \code{NULL}, in which case such 
#'structures will not be taken into account. 
#'
#'@param Sigma_xi a matrix of size \code{K x K} containing the covariance matrix
#'of the \code{K} random effects corresponding to \code{phi}.
#'
#'@param na.rm logical: should missing values (including \code{NA} and
#'\code{NaN}) be omitted from the calculations? Default is \code{FALSE}.
#'
#'@return A list with the following elements:\itemize{
#'   \item \code{score}: approximation of the set observed score
#'   \item \code{q}: observation-level contributions to the score
#'   \item \code{q_ext}: "pseudo-observations" used to compute the covariance of q,
#'    taking into account the contributions of OLS estimates
#'   \item \code{gene_scores_unscaled}: a vector of the approximations of the
#'   individual gene scores
#' }
#'
#'@examples
#'set.seed(123)
#'
#'##generate some fake data
#'########################
#'n <- 100
#'r <- 12
#'t <- matrix(rep(1:3), r/3, ncol=1, nrow=r)
#'sigma <- 0.4
#'b0 <- 1
#'
#'#under the null:
#'b1 <- 0
#'#under the alternative:
#'b1 <- 0.7
#'y.tilde <- b0 + b1*t + rnorm(r, sd = sigma)
#'y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
#'       matrix(rep(y.tilde, n), ncol=n, nrow=r))
#'x <- matrix(1, ncol=1, nrow=r)
#'
#'#run test
#'scoreTest <- vc_score(y, x, phi=t, w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
#'                     Sigma_xi=matrix(1), indiv=rep(1:(r/3), each=3))
#'scoreTest$score
#'
#'@importFrom stats model.matrix
#'@importFrom reshape2 melt
#'@importFrom Matrix bdiag
#'@export
vc_score <- function(y, x, indiv = c(1:ncol(y)), phi, w = NULL, Sigma_xi = NULL,
                     Sigma = NULL, na.rm = FALSE) {
  
  ## dimensions, formatting and validity checks---
  
  if (!is.null(Sigma)){
    cor_structure = TRUE # set argument as TRUE if Sigma provided
  } else {
    cor_structure = FALSE
  }
  
  # Check finiteness of covariance matrix entries
  if (cor_structure == TRUE & (sum(!is.finite(unlist(Sigma)))) > 0) {
    stop("At least 1 non-finite covariance entry in Sigma")
  }
  # Check non-singularity of covariance matrices
  if(cor_structure == TRUE & any(unlist(lapply(Sigma, 
                                               FUN = function(x){any(svd(x)$d < 1e-10)})))){
    cor_structure = FALSE #If any Sigma matrices are singular, we do not use them
  }
  
  # Check data is provided in matrix format
  stopifnot(is.matrix(y)) 
  stopifnot(is.matrix(x))
  stopifnot(is.matrix(phi))
  
  # Set arguments based on input dimensions
  p <- nrow(y)  # the number of genes in the set
  n <- ncol(y)  # the number of samples measured
  q <- ncol(x)  # the number of covariates
  K <- ncol(phi) # the number of test variables
  
  # converting indiv vector to numeric factor
  indiv <- as.factor(as.numeric(as.factor(indiv))) 
  n_indiv <- length(levels(indiv)) # number of individuals
  n_i = (summary(indiv, maxsum = n_indiv) %>% as.vector()) # number of observations per individual
  
  # Check concordance of input matrices
  stopifnot(nrow(x) == n) # covariate matrix must have samples as rows
  stopifnot(nrow(phi) == n) # test variable matrix must have samples as rows
  stopifnot(length(indiv) == n) # vector of sample indices must have n elements
  
  # Check heteroskedasticity weights
  if (is.null(w)){
    w = matrix(1, nrow = p, ncol = n) # if no weights specified, heteroskedasticity not accounted for
  } else {
    stopifnot(nrow(w) == p | ncol(w) == n) # else weights must have column samples and row genes
  }
  
  # Check finiteness of weights
  if (sum(!is.finite(w)) > 0) {
    stop("At least 1 non-finite weight in 'w'") 
  }
  
  # check correct number and dimensions of covariance matrices
  if (cor_structure == TRUE){
    stopifnot(length(Sigma) == n_indiv) 
    stopifnot(all(unlist(lapply(lapply(Sigma, dim), `[[`, 1)) == n_i*p))
  }
  
  # Check random effects covariance matrix
  if (is.null(Sigma_xi)){
    Sigma_xi = diag(ncol(phi)) # If not provided, set as identity
  }
  stopifnot(ncol(Sigma_xi) == K) # Check correct dimensions
  
  ## OLS for conditional mean -----
  if (na.rm & sum(is.na(y)) > 0) {
    y_0 <- y
    y_0[is.na(y_0)] <- 0 } else { # removing NAs if argument selected
      y_0 <- y
    }
  
  y_T = t(y_0) # Transpose of expression matrix
  # OLS estimate of mean effect, assuming homogeneous effect of covariates across genes
  alpha <- solve(crossprod(x)) %*% t(x) %*% rowMeans(y_T, na.rm = na.rm)
  # Centered gene expression given OLS estimate
  yt_mu <- y_T - do.call(cbind, replicate(p, x %*% alpha, simplify = FALSE))
  y_mu = t(yt_mu) # Transpose of centered gene expression
  
  if (cor_structure == TRUE){ # Case for when correlation structures are estimated
    
    # List of phi matrices in block diagonal format
    # First create list of values split by individual index, then apply bdiag function with lapply
    phi_diag_list = lapply(split(phi, indiv), FUN = function(input_list){
      as.matrix(bdiag(rep(list(input_list), p)))
    })
    
    # Transforming to long format - i.e. each column an individual with "stacked" observations
    y_mu_long = melt(cbind(data.frame(yt_mu), individual = indiv), id.vars = "individual")
    
    # Creating a list with ith element as the (long format) centered expression for individual i
    y_mu_long_list = split(y_mu_long$value, y_mu_long$individual)
    
    # Inverting covariance matrices using Cholesky decomposition (faster than solve())
    sigma_inv_list = lapply(lapply(Sigma, chol), Matrix::chol2inv)
    
    # mapply to perform the required matrix multiplication (faster than for loop)
    Q_indiv = mapply(FUN = function(A,B,C){return(crossprod(A,B) %*% C)}, 
                     y_mu_long_list, sigma_inv_list, phi_diag_list)
    
    # Final outputs 
    q = t(Q_indiv) # Individual level test statistic contributions 
    qq <- colSums(q, na.rm = na.rm)^2/n_indiv  # Genewise, variable specific scores
    gene_Q <- rowSums(matrix(qq, ncol = K))  # Gene scores
    Q <- sum(qq) # set score
    
    ## Calculating "pseudo-observations" for estimating the covariance of q
    sig_xi_sqrt <- (Sigma_xi * diag(K)) ^ (-0.5)
    sig_eps_inv_T <- t(w)
    
    phi_sig_xi_sqrt <- phi %*% sig_xi_sqrt # =phi
    
    T_fast <- do.call(cbind, replicate(K, sig_eps_inv_T, simplify = FALSE)) *
      matrix(apply(phi_sig_xi_sqrt, 2, rep, p), ncol = p * K)
    
    XT_fast <- crossprod(x, T_fast)/n_indiv  # Covariate contributions?
    avg_xtx_inv_tx <- n_indiv * tcrossprod(solve(crossprod(x, x)), x)
    U_XT <- matrix(yt_mu, ncol = p * K, nrow = n) *
      crossprod(avg_xtx_inv_tx, XT_fast)
    if (na.rm & sum(is.na(U_XT)) > 0) {
      U_XT[is.na(U_XT)] <- 0
    }
    
    if (length(levels(indiv)) > 1) {
      indiv_mat <- stats::model.matrix(~0 + factor(indiv))
    } else {
      indiv_mat <- matrix(as.numeric(indiv), ncol = 1)
    }
    
    U_XT_indiv <- crossprod(indiv_mat, U_XT)
    q_ext <- q - U_XT_indiv #pseudo_observations
    
    return(list(score = Q, q = q, q_ext = q_ext, gene_scores_unscaled = gene_Q))
    
  } else { # Case for when correlation structures are NOT estimated
    sig_xi_sqrt <- (Sigma_xi * diag(K)) %^% (-0.5)
    sig_eps_inv_T <- t(w)
    phi_sig_xi_sqrt <- phi %*% sig_xi_sqrt
    
    T_fast <- do.call(cbind, replicate(K, sig_eps_inv_T, simplify = FALSE)) *
      matrix(apply(phi_sig_xi_sqrt, 2, rep, p), ncol = p * K)
    ##---------------------
    ## the structure of T_fast is time_basis_1*gene_1, time_basis_1*gene_2, ...,
    ## time_basis_1*gene_p, ..., time_basis_K*gene_1, ..., time_basis_K*gene_p
    ##----------------------------
    q_fast <- matrix(yt_mu, ncol = p * K, nrow = n) * T_fast
    if (length(levels(indiv)) > 1) {
      indiv_mat <- stats::model.matrix(~0 + factor(indiv))
    } else {
      indiv_mat <- matrix(as.numeric(indiv), ncol = 1)
    }
    
    if (na.rm & sum(is.na(q_fast)) > 0) {
      q_fast[is.na(q_fast)] <- 0
    }
    q <- crossprod(indiv_mat, q_fast)
    XT_fast <- crossprod(x, T_fast)/n_indiv
    avg_xtx_inv_tx <- n_indiv * tcrossprod(solve(crossprod(x, x)), x)
    U_XT <- matrix(yt_mu, ncol = p * K, nrow = n) *
      crossprod(avg_xtx_inv_tx, XT_fast)
    if (na.rm & sum(is.na(U_XT)) > 0) {
      U_XT[is.na(U_XT)] <- 0
    }
    
    U_XT_indiv <- crossprod(indiv_mat, U_XT)
    q_ext <- q - U_XT_indiv
    qq <- colSums(q, na.rm = na.rm)^2/n_indiv  # genewise, variable specific scores
    gene_Q <- rowSums(matrix(qq, ncol = K))  # gene scores
    
    Q <- sum(qq) # set score
    
    return(list(score = Q, q = q, q_ext = q_ext, gene_scores_unscaled = gene_Q))
  } 
}
