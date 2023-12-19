#'Asymptotic variance component test statistic and p-value
#'
#'This function computes an approximation of the variance component test based
#'on the asymptotic distribution of a mixture of \eqn{\chi^{2}}s using the saddlepoint
#'method from \code{\link[survey]{pchisqsum}}, as per Chen & Lumley 20219 CSDA.
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
#'@param genewise_pvals a logical flag indicating whether gene-wise p-values
#'should be returned. Default is \code{FALSE} in which case gene set p-value is
#'computed and returned instead.
#'
#'@param homogen_traj a logical flag indicating whether trajectories should be
#'considered homogeneous. Default is \code{FALSE} in which case trajectories
#'are not only tested for trend, but also for heterogeneity.
#'
#'@param na.rm logical: should missing values (including \code{NA} and
#'\code{NaN}) be omitted from the calculations? Default is \code{FALSE}.
#'
#'@return A list with the following elements when the set p-value is computed:
#'\itemize{
#'   \item \code{set_score_obs}: the approximation of the observed set score
#'   \item \code{set_pval}: the associated set p-value
#' }
#' or a list with the following elements when gene-wise p-values are computed:
#' \itemize{
#'   \item \code{gene_scores_obs}: vector of approximating the observed
#'   gene-wise scores
#'   \item \code{gene_pvals}: vector of associated gene-wise p-values
#' }
#'
#'
#'@seealso \code{\link[survey]{pchisqsum}}
#'
#'@references
#'Chen T & Lumley T (2019), Numerical evaluation of methods approximating the
#'distribution of a large quadratic form in normal variables, Computational
#'Statistics & Data Analysis, 139:75-81.
#'
#'@examples
#'set.seed(123)
#'
#'##generate some fake data
#'########################
# n <- 100
# r <- 12
# t <- matrix(rep(1:(r/4)), 4, ncol=1, nrow=r)
# sigma <- 0.4
# b0 <- 1
# 
# #under the null:
# b1 <- 0
# #under the alternative:
# #b1 <- 0.5
# y.tilde <- b0 + b1*t + rnorm(r, sd = sigma)
# y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
#       matrix(rep(y.tilde, n), ncol=n, nrow=r))
# x <- matrix(1, ncol=1, nrow=r)
#'
#'#run test
#'asymTestRes <- vc_test_asym(y, x, phi=cbind(t, t^2),
#'                            w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
#'                            Sigma_xi=diag(2), indiv=1:r, genewise_pvals=TRUE)
#'plot(density(asymTestRes$gene_pvals))
#'quantile(asymTestRes$gene_pvals)
#'
#'@importFrom survey pchisqsum
#'@importFrom CompQuadForm davies
#'@importFrom stats pchisq cov
#'@importFrom matrixStats colVars
#'
#'@export
vc_test_asym <- function(y, x, indiv = rep(1, nrow(x)), phi, w = NULL,
                         Sigma_xi = diag(ncol(phi)),
                         Sigma = NULL, 
                         genewise_pvals = FALSE, homogen_traj = FALSE,
                         na.rm = FALSE) {
    ## dimensions, formatting and validity checks------
    
    # check data is provided in matrix format
    stopifnot(is.matrix(y)) 
    stopifnot(is.matrix(x))
    stopifnot(is.matrix(phi))
    
    # Logical flag if sigma matrices given
    if (!is.null(Sigma)){
      cor_structure = TRUE 
    } else {
      cor_structure = FALSE
    }
  
    #Checking for finiteness of Sigma matrix entries 
    if (cor_structure == TRUE & (sum(!is.finite(unlist(Sigma)))) > 0) {
      stop("At least 1 non-finite covariance entry in Sigma")
    }

    #Checking for singularity of Sigma matrices using SVD
    if(cor_structure == TRUE & any(unlist(lapply(Sigma, 
                                                FUN = function(x){any(svd(x)$d < 1e-10)})))){
      message("At least one singular covariance matrix for a geneset. \n",
              "Will not incorporate correlation structures into estimation for this gene set. \n")
      cor_structure = FALSE 
    }
  
    p <- nrow(y)  # the number of genes in the set
    n <- ncol(y)  # the number of samples measured
    q <- ncol(x)  # the number of covariates
    K <- ncol(phi) # the number of test variables
  
    indiv <- as.factor(as.numeric(as.factor(indiv))) # converting indiv vector to numeric factor
    n_indiv <- length(levels(indiv)) # number of individuals
    n_i = (summary(indiv ,maxsum = n_indiv) %>% as.vector()) # number of observations per individual
  
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
    
    # Check correct number and dimensions of covariance matrices 
    if (cor_structure == TRUE){
      stopifnot(length(Sigma) == n_indiv)
      stopifnot(all(unlist(lapply(lapply(Sigma, dim), `[[`, 1)) == n_i*p))
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
    
    # Call score function
    if (homogen_traj) {
        score_list <- vc_score_h(y = y, x = x, indiv = factor(indiv), phi = phi,
                                 w = w, Sigma_xi = Sigma_xi, na.rm = na.rm)
    } else {
        score_list <- vc_score(y = y, x = x, indiv = factor(indiv), phi = phi,
                               w = w, Sigma = Sigma , Sigma_xi = Sigma_xi, na.rm = na.rm)
    }

    if (p * n_indiv < 1) {
        stop("no gene measured/no sample included ...")
    }
    
    Q_indiv = t(score_list$q)
    
    if (genewise_pvals) {
        gene_scores_obs <- score_list$gene_scores_unscaled
        if (n_indiv == 1) {
            pv <- stats::pchisq(gene_scores_obs, df = 1, lower.tail = FALSE)
        } else if (K == 1) {
          gene_lambda <- matrixStats::colVars(score_list$q_ext)
          #gene_lambda <- matrixStats::colVars(score_list$q)
          pv <- stats::pchisq(gene_scores_obs/gene_lambda, df = 1,
                              lower.tail = FALSE)
        } else {
            gene_inds <- lapply(seq_len(p), function(x) {
                x + (p) * (seq_len(K) - 1)
            })

            gene_lambda <- lapply(gene_inds, function(x) {
                Sig_q_gene <- cov(Q_indiv[, x, drop = FALSE])
                lam <- tryCatch(eigen(Sig_q_gene, symmetric = TRUE, only.values = TRUE)$values,
                                error=function(cond){return(NULL)}
                )
                if (is.null(lam)){
                    lam <- tryCatch(svd(Sig_q_gene)$d,
                                    error=function(cond){return(NULL)}
                    )
                }
                if (is.null(lam)){
                    lam <- tryCatch(svd(Sig_q_gene)$d,
                                    error=function(cond){
                                        warning("SVD decomposition failed for at least one ",
                                                "gene")
                                        return(NA)
                                    })
                }
                return(lam)
            })
            pv <- try(mapply(FUN = survey::pchisqsum,
                             x = gene_scores_obs,
                             df=1,
                             a = gene_lambda,
                             lower.tail=FALSE,
                             method = "saddlepoint"
            ), silent = TRUE)
            if(inherits(pv, "try-error")){
              ## old slow CompQuadForm::davies method  (sometimes accuracy issues with low p-vals)
              pv <- unlist(mapply(FUN = CompQuadForm::davies,
                                  q = gene_scores_obs,
                                  lambda = gene_lambda, lim = 15000,
                                  acc = 5e-04)["Qq", ])
            }
        }

        names(pv) <- rownames(y)

        ans <- list("gene_scores_obs" = gene_scores_obs,
                    "gene_pvals" = pv)

    } else {

        if (n_indiv == 1) {
            Gamma <- matrix(1, p, p)
        } else {
          # Q_bar = (1/n_indiv)*rowSums(Q_indiv) # mean vector 
          # Q_list = vector("list", n_indiv) #list to store p times p individual contritubion matrices
          # for (i in c(1:n_indiv)){
          #   Q_list[[i]] = tcrossprod((Q_indiv[,i] - Q_bar))
          # }
          # Gamma = (1/n_indiv)*Reduce('+', Q_list) # scaled estimated covariance matrix of 
          #                                         # individual-level contributions
          
          Gamma = cov(score_list$q_ext)
        }

        lam <- tryCatch(eigen(Gamma, symmetric = TRUE, only.values = TRUE)$values,
                        error=function(cond){return(NULL)} #estimated eigenvalues
        )
        if (is.null(lam)) {
            lam <- tryCatch(svd(Gamma)$d,
                            error=function(cond){return(NULL)}
            )
        }
        if (is.null(lam)) {
            lam <- tryCatch(svd(round(Gamma, 6))$d,
                            error=function(cond){
                                stop("SVD decomposition failed")
                            })
        }

        dv <- survey::pchisqsum(x = score_list$score,
                                    df=1,
                                    a = lam,
                                    lower.tail = FALSE,
                                    method = "saddlepoint") # estimated p value

        ans <- list("set_score_obs" = score_list$score,
                    "set_pval" = dv)
    }

    return(ans)
}

