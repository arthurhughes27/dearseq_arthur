#'Computes variance component score test statistics
#'
#'This function computes the variance component score test statistics
#'
#'@keywords internal
#'
#'@param y a numeric matrix of dim \code{p x n} containing the raw RNA-seq
#'counts for \code{p} genes from \code{n} samples.
#'
#'@param x a numeric design matrix of dim \code{n x q} containing the \code{q}
#'covariates to be adjusted for.
#'
#'@param indiv a vector of length \code{n} containing the information for
#'attributing each sample to one of the studied individuals. Coerced
#'to be a \code{factor}. Default is \code{NULL}, in which case it is assumed that
#'each row represents a distinct individual.
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the
#'\code{K} variables to be tested
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
#'of the \code{K} random effects on \code{phi}.
#'
#'@param na_rm logical: should missing values (including \code{NA} and
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
#'
#'@export
vc_score <- function(y, x, indiv = NULL, phi, w = NULL, Sigma_xi = diag(ncol(phi)),
                     Sigma = NULL, na_rm = FALSE) {
  
    ## dimensions, formatting and validity checks------
    if (sum(!is.finite(w)) > 0) {
        stop("At least 1 non-finite weight in 'w'")
    }
  
    if (!is.null(Sigma)){
      cor_structure = TRUE # set argument if Sigma provided
    }
  
    if (cor_structure == TRUE & (sum(!is.finite(unlist(Sigma)))) > 0) {
      stop("At least 1 non-finite covariance entry in Sigma")
    }

    stopifnot(is.matrix(y))
    stopifnot(is.matrix(x))
    stopifnot(is.matrix(phi))
    
    p <- nrow(y)  # the number of genes in the set
    n <- ncol(y)  # the number of samples measured
    q <- ncol(x)  # the number of covariates
    K <- ncol(phi) # the number of test variables
    
    if (is.null(indiv)){
      indiv = c(1:n) # set argument if Sigma provided
    }
    
    indiv <- as.factor(as.numeric(as.factor(indiv))) #converting indiv vector to numeric factor
    n_indiv <- length(levels(indiv)) #number of individuals
    n_i = (summary(indiv) %>% as.vector()) #number of observations per individual
    
    stopifnot(nrow(x) == n) # x must have samples as rows
    stopifnot(length(indiv) == n) # vector of sample indices must have n elements
    
    if (is.null(w)){
      w = matrix(1, nrow = p, ncol = n) # if no weights specified, heteroskedasticity not accounted for
    } else {
      stopifnot(nrow(w) == p | ncol(w) == n) #else weights must have column samples and row genes
    }
    stopifnot(nrow(phi) == n) # test variable must have samples as rows
    stopifnot(length(indiv) == n) # vector of sample indices must have n elements
    if (cor_structure == TRUE){
      stopifnot(length(Sigma) == n_indiv) # check correct number of covariance matrices
      stopifnot(all(unlist(lapply(lapply(Sigma, dim), `[[`, 1)) == n_i*p))
      # check correct dimensions of all covariance matrices
    }
    
    # the number of random effects
    if (length(Sigma_xi) == 1) {
        K <- 1
        Sigma_xi <- matrix(Sigma_xi, K, K)
    } else {
        K <- nrow(Sigma_xi)
        stopifnot(ncol(Sigma_xi) == K)
    }
    stopifnot(K == K)

    ## Generalised OLS for conditional mean -----
    if (na_rm & sum(is.na(y)) > 0) {
      y_0 <- y
      y_0[is.na(y_0)] <- 0 } else {
      y_0 <- y
      }
    
    sigma_inv = vector("list", n_indiv)
    y_mu_list = vector("list", n_indiv)
    for (i in c(1:n_indiv)){
      if (n_i[i] == 1){
        X_diag = t(as.matrix(bdiag(rep(list(x[indiv == i,]), p))))
      } else {
        X_diag = as.matrix(bdiag(rep(list(x[indiv == i,]), p)))
      }
    
      yi = rep(0,n_i[i]*p)
      for (g in c(1:p)){
        yi[c(((g-1)*n_i[i]+1):(g*n_i[i]))] = y_0[g,c(indiv == i)]
      }
    
      if (cor_structure == TRUE){
        sigma_inv[[i]] = solve(sigma_list[[i]]) # full unstructured covariance matrix
      } else {
        wcol = rep(0,n_i[i]*p)
        for (g in c(1:p)){
          wcol[c(((g-1)*n_i[i]+1):(g*n_i[i]))] = w[g,indiv == i]
        }
        sigma_inv[[i]] = diag(wcol)
      }
    
      y_mu_list[[i]] = rep(0, n_i[i]*p)
      xt_sig_x_inv = Ginv((t(X_diag) %*% sigma_inv[[i]] %*% X_diag))
      y_mu_list[[i]] = yi - X_diag %*% xt_sig_x_inv %*% t(X_diag) %*% sigma_inv[[i]] %*% yi
      }
    
    
    Q_mat = matrix(0,nrow = p*K, ncol = n_indiv)
    
    for (i in c(1:n_indiv)){
    if (n_i[i] == 1){
      phi_diag = t(as.matrix(bdiag(rep(list(phi[indiv == i,]), p))))
    } else {
      phi_diag = as.matrix(bdiag(rep(list(phi[indiv == i,]), p)))
    }
    
    Q_mat[,i] = t(y_mu_list[[i]]) %*% sigma_inv[[i]] %*% phi_diag}
    
    Q_bar = (1/n_indiv)*rowSums(Q_mat)
    Q_list = vector("list", n_indiv)
    for (i in c(1:n_indiv)){
      Q_list[[i]] = tcrossprod(Q_mat[,i] - Q_bar)
    }
    Gamma = Reduce('+', Q_list)
    
    q_col = n^{-1/2}*rowSums(Q_mat)
    Q = crossprod(q_col)
    gene_Q = rowSums(Q_mat) # pay attention to case with >1 testing variable 
    
    return(list(score = Q, q = Q_mat, q_ext = Q_mat,
                gene_scores_unscaled = gene_Q, Gamma = Gamma))
    
    # Blocked out by arthur
    # sig_xi_sqrt <- (Sigma_xi * diag(K)) ^ (-0.5)
    # sig_eps_inv_T <- t(w)
    # 
    # phi_sig_xi_sqrt <- phi %*% sig_xi_sqrt
    # 
    # 
    # 
    # T_fast <- do.call(cbind, replicate(K, sig_eps_inv_T, simplify = FALSE)) *
    #     matrix(apply(phi_sig_xi_sqrt, 2, rep, p), ncol = p * K)
    # ##---------------------
    # ## the structure of T_fast is time_basis_1*gene_1, time_basis_1*gene_2, ...,
    # ## time_basis_1*gene_p, ..., time_basis_K*gene_1, ..., time_basis_K*gene_p
    # ##----------------------------
    # q_fast <- matrix(yt_mu, ncol = p * K, nrow = n) * T_fast
    # 
    # ## dplyr seems to be less efficient here
    # ## q_fast_tb <- tibble::as_tibble(cbind.data.frame(indiv, q_fast))
    # ## q_dp <- q_fast_tb %>% group_by(indiv) %>% summarise_all(sum)
    # 
    # ## aggregate is much slower also
    # ## qtemp <- aggregate(. ~ indiv, cbind.data.frame(indiv, q_fast), sum)
    # ## qtemp <- aggregate(. ~ indiv, cbind.data.frame(indiv, q_fast), sum)
    # 
    # ## data.table hard to test, but seems to be at least 10 times slower on big
    # ## datasets (weird)
    # ## m_dt <- data.table('indiv'=factor(rep(c(1:20), each=5)), mbig)
    # ## temp <- m_dt[, lapply(.SD, sum), by=indiv]
    # 
    # 
    # # the 2 'by' statements below used to represent the longest AND most memory
    # # intensive part of this for genewise analysis:
    # if (length(levels(indiv)) > 1) {
    #     indiv_mat <- stats::model.matrix(~0 + factor(indiv))
    # } else {
    #     indiv_mat <- matrix(as.numeric(indiv), ncol = 1)
    # }
    # 
    # if (na_rm & sum(is.na(q_fast)) > 0) {
    #     q_fast[is.na(q_fast)] <- 0
    # }
    # q <- crossprod(indiv_mat, q_fast)
    # XT_fast <- crossprod(x, T_fast)/n_indiv
    # avg_xtx_inv_tx <- n_indiv * tcrossprod(solve(crossprod(x, x)), x)
    # U_XT <- matrix(yt_mu, ncol = p * K, nrow = n) *
    #     crossprod(avg_xtx_inv_tx, XT_fast)
    # if (na_rm & sum(is.na(U_XT)) > 0) {
    #     U_XT[is.na(U_XT)] <- 0
    # }
    # U_XT_indiv <- crossprod(indiv_mat, U_XT)
    # q_ext <- q - U_XT_indiv 
    # # sapply(1:6, function(i){(q_ext[i,] - q_ext_fast_indiv[i,])})
    # 
    # 
    # 
    # qq <- colSums(q, na.rm = na_rm)^2/n_indiv  # genewise scores
    # 
    # ##qq <- unlist(by(data=matrix(qq, ncol=1), INDICES=rep(1:p, K), FUN=sum,
    # ##                simplify = FALSE)) # veryslow
    # ## gene_inds <- lapply(1:p, function(x){x + (p)*(0:(K-1))})
    # ## gene_Q <- sapply(gene_inds, function(x) sum(qq[x])) # old computation
    # ## gene_Q <- tcrossprod(qq, matrix(diag(p), nrow=p, ncol=p*K))[1, ] # faster
    # gene_Q <- rowSums(matrix(qq, ncol = K))  # even faster
    # 
    # QQ <- sum(qq)  #n_indiv=nrow(q) # set score

    # return(list(score = QQ, q = q, q_ext = q_ext,
    #             gene_scores_unscaled = gene_Q))
}
