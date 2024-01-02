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
#'   \item \code{q_ext}: pseudo-observations used to compute the covariance,
#'    taking into account the contributions of OLS estimates
#'   \item \code{gene_scores_unscaled}: a vector of the approximations of the
#'   individual gene scores
#' }
#'
#'
#'@examples
#'set.seed(123)
#'
#'##generate some fake data
#'########################
#'ng <- 100
#'nindiv <- 30
#'nt <- 5
#'nsample <- nindiv*nt
#'tim <- matrix(rep(1:nt), nindiv, ncol=1, nrow=nsample)
#'tim2 <- tim^2
#'sigma <- 5
#'b0 <- 10
#'
#'#under the null:
#'beta1 <- rnorm(n=ng, 0, sd=0)
#'#under the (heterogen) alternative:
#'beta1 <- rnorm(n=ng, 0, sd=0.1)
#'#under the (homogen) alternative:
#'beta1 <- rnorm(n=ng, 0.06, sd=0)
#'
#'y.tilde <- b0 + rnorm(ng, sd = sigma)
#'y <- t(matrix(rep(y.tilde, nsample), ncol=ng, nrow=nsample, byrow=TRUE) +
#'       matrix(rep(beta1, each=nsample), ncol=ng, nrow=nsample, byrow=FALSE) *
#'            matrix(rep(tim, ng), ncol=ng, nrow=nsample, byrow=FALSE) +
#'       #matrix(rep(beta1, each=nsample), ncol=ng, nrow=nsample, byrow=FALSE) *
#'       #    matrix(rep(tim2, ng), ncol=ng, nrow=nsample, byrow=FALSE) +
#'       matrix(rnorm(ng*nsample, sd = sigma), ncol=ng, nrow=nsample,
#'              byrow=FALSE)
#'       )
#'myindiv <- rep(1:nindiv, each=nt)
#'x <- cbind(1, myindiv/2==floor(myindiv/2))
#'myw <- matrix(rnorm(nsample*ng, sd=0.1), ncol=nsample, nrow=ng)
#'
#'#run test
#'score_homogen <- vc_score_h(y, x, phi=tim, indiv=myindiv,
#'                            w=myw, Sigma_xi=cov(tim))
#'score_homogen$score
#'
#'score_heterogen <- vc_score(y, x, phi=tim, indiv=myindiv,
#'                            w=myw, Sigma_xi=cov(tim))
#'score_heterogen$score
#'
#'scoreTest_homogen <- vc_test_asym(y, x, phi=tim, indiv=rep(1:nindiv, each=nt),
#'                                  w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
#'                                  Sigma_xi=cov(tim),
#'                                  homogen_traj = TRUE)
#'scoreTest_homogen$set_pval
#'scoreTest_heterogen <- vc_test_asym(y, x, phi=tim, indiv=rep(1:nindiv,
#'                                                          each=nt),
#'                                    w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
#'                                    Sigma_xi=cov(tim),
#'                                    homogen_traj = FALSE)
#'scoreTest_heterogen$set_pval
#'@importFrom stats model.matrix
#'@importFrom reshape2 melt
#'@importFrom Matrix bdiag
#'@export
vc_score_h <- function(y, x, indiv = c(1:ncol(y)), phi, w = NULL, Sigma_xi = diag(ncol(phi)),
                     Sigma = NULL, na.rm = FALSE) {

    ## dimensions, formatting and validity checks---
  
    if (!is.null(Sigma)){
     cor_structure = TRUE # set argument if Sigma provided
    } else {
     cor_structure = FALSE
    }
  
    if (cor_structure == TRUE & (sum(!is.finite(unlist(Sigma)))) > 0) {
     stop("At least 1 non-finite covariance entry in Sigma") # check finite covariance entries
    }
  
    if(cor_structure == TRUE & any(unlist(lapply(Sigma, 
                                                 FUN = function(x){any(svd(x)$d < 1e-10)})))){
     cor_structure = FALSE #If any Sigma matrices are singular, do not use them
    }
  
  
    stopifnot(is.matrix(y)) # check data is provided in matrix format
    stopifnot(is.matrix(x))
    stopifnot(is.matrix(phi))
  
    p <- nrow(y)  # the number of genes in the set
    n <- ncol(y)  # the number of samples measured
    q <- ncol(x)  # the number of covariates
    K <- ncol(phi) # the number of test variables
  
    indiv <- as.factor(as.numeric(as.factor(indiv))) # converting indiv vector to numeric factor
    n_indiv <- length(levels(indiv)) # number of individuals
    n_i = (summary(indiv) %>% as.vector()) # number of observations per individual
  
    stopifnot(nrow(x) == n) # covariate matrix must have samples as rows
    stopifnot(nrow(phi) == n) # test variable matrix must have samples as rows
    stopifnot(length(indiv) == n) # vector of sample indices must have n elements
  
    if (is.null(w)){
      w = matrix(1, nrow = p, ncol = n) # if no weights specified, heteroskedasticity not accounted for
    } else {
     stopifnot(nrow(w) == p | ncol(w) == n) # else weights must have column samples and row genes
    } 
  
  # Check finiteness of weights
    if (sum(!is.finite(w)) > 0) {
      stop("At least 1 non-finite weight in 'w'") 
    }
  
    if (cor_structure == TRUE){
      stopifnot(length(Sigma) == n_indiv) # check correct number of covariance matrices
     stopifnot(all(unlist(lapply(lapply(Sigma, dim), `[[`, 1)) == n_i*p))
     # check correct dimensions of all covariance matrices
    }
  
    # the number of random effects
    if (is.null(Sigma_xi)){
     Sigma_xi = diag(ncol(phi))
    }
    if (length(Sigma_xi) == 1) {
     K <- 1
     Sigma_xi <- matrix(Sigma_xi, K, K)
    } else {
     K <- nrow(Sigma_xi)
     stopifnot(ncol(Sigma_xi) == K)
    }
    stopifnot(K == K)
  
    ## OLS for conditional mean -----
    if (na.rm & sum(is.na(y)) > 0) {
     y_0 <- y
      y_0[is.na(y_0)] <- 0 } else { # removing NAs if argument selected
        y_0 <- y
      }
    
    y_T = t(y_0)
    alpha <- solve(crossprod(x)) %*% t(x) %*% rowMeans(y_T, na.rm = na.rm)
    yt_mu <- y_T - do.call(cbind, replicate(p, x %*% alpha, simplify = FALSE))
    y_mu = t(yt_mu)

    ## test statistic computation ------
    sig_xi_sqrt <- (Sigma_xi * diag(K))^(-0.5)
    
    sigma_inv = vector("list", n_indiv) # List to store the inverted covariance matrices
    Q_indiv = matrix(0,nrow = p*K, ncol = n_indiv) #matrix for individual level contributions
    for (i in c(1:n_indiv)){ 
      # Test variable matrix in block diag form
      if (n_i[i] == 1){
        phi_diag = t(as.matrix(bdiag(rep(list(phi[indiv == i,]), p))))
      } else {
        phi_diag = as.matrix(bdiag(rep(list(phi[indiv == i,]), p)))
      }
      
      y_mui = rep(0,n_i[i]*p) # temporary column vector for centered gene expression
      for (g in c(1:p)){
        y_mui[c(((g-1)*n_i[i]+1):(g*n_i[i]))] = y_mu[g,c(indiv == i)] #rearranging samples 
                                                                      # to long format
      } 
      
      if (cor_structure == TRUE){
        sigma_inv[[i]] = solve(Sigma[[i]]) # inverse of full unstructured covariance matrix if given
      } else {
        wcol = rep(0,n_i[i]*p) #temporary column vector to store obs-level weights
        for (g in c(1:p)){
          wcol[c(((g-1)*n_i[i]+1):(g*n_i[i]))] = w[g,indiv == i] # long format 
        }
        sigma_inv[[i]] = wcol %>% as.vector() %>% diag() 
        # covariance matrix is simply weights if no correlation structure desired
        # (weights are already inverse variances)
      }
      
      Q_indiv[,i] = crossprod((y_mui),sigma_inv[[i]]) %*% phi_diag 
    }
    
    q_col = n_indiv^{-1/2}*rowSums(Q_indiv)
    Q = crossprod(q_col)
    gene_Q_unscaled = rowSums(Q_indiv) # pay attention to case with >1 testing variable 
    
    return(list(score = Q, Q_indiv = Q_indiv, gene_scores_unscaled = gene_Q_unscaled))
}
