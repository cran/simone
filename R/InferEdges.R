## __________________________________________________________
##
## InferEdges
##
## INPUT
##	data	: n*p empirical data matrix (e.g. gene expressions). 
##		  NA values are ignored.
##	penalty : penalty mask to use (penalty.handler object,
##		  numeric matrix, numeric vector or scalar)
##      method : a string that contains the inference method to
##                perform, either "glasso", "regressionAND" or
##                "regressionOR" (the latters correspond to
##                Meinshausen and Buhlman's approximation, with AND or
##                OR rule to symmetrize the result). Default is
##                "glasso".
##
## OPTIONAL INPUT
##	eps 	: threshold for covergence
##	maxIt 	: maximum number of iterations for outer loop
##
## OUTPUT
##     K.hat        : estimate of the inverse of the covariance matrix
##     Sigma.hat    : estimate of the covariance matrix
##
## Estimate the inverse covariance matrix from the empirical
## covariance matrix solving p independent l1-penalized regressions
## (Meinshausen and Bulhman) or the a l1-penalized likelihood
## problem (Banerjee et al, Friedman et al).
## __________________________________________________________
##
InferEdges <- function(data, penalty=2/(nrow(data)*BaseLambdaValue(data,risk=0.3)),
                            method = "glasso", ...) {

  ## defaults for hidden optional parameters ##
  eps 	<- sub.param("eps"  , 1e-12, ...)
  maxIt <- sub.param("maxIt", 1e4, ...)
  Sigma.hat <- sub.param("Sigma.hat", NULL, ...)

  ## Empirical covariance matrix ##
  S <- var(data, na.rm=TRUE)
  p <- nrow(S)

  if (is.matrix(penalty)) {
    Rho <- penalty
    Rho[is.infinite(Rho)] <- -1
  } else  if (length(penalty)==1) {
    Rho <- matrix(penalty, p, p)
  } else {
    return(NULL)
  }
    
  ## Initializing the Matrice ##
  K.hat <- matrix(0, nrow=p,ncol=p)
  
  ## SOLVING PENALIZED PSEUDO-LIKELIHOOD, SYMMETRIZE WITH AND RULE ##
  if ("regressionAND" %in% method) {
    out <- .C("regLasso",
              as.integer (p),
              as.integer(0), # AND rule: boolean fixed to 0 
              as.double  (S),
              as.double  (Rho),
              as.double  (eps),
              K.hat =  as.double(K.hat),
              PACKAGE="simone")
    
    K.hat <- matrix(out$K.hat, ncol=p)
    # Re-normalize the solution to get the associated precision matrix #
    K.hat <- Symmetrize(- diag(1/diag(S)) %*% (K.hat * upper.tri(K.hat))
                        + diag(1/diag(S)) )
  }
  
  ## SOLVING PENALIZED PSEUDO-LIKELIHOOD, SYMMETRIZE WITH OR RULE ##
  else if  ("regressionOR" %in% method) {
    
    out <- .C("regLasso",
              as.integer (p),
              as.integer(1), # OR rule: boolean fixed to anything but 0 
              as.double  (S),
              as.double  (Rho),
              as.double  (eps),
              K.hat =  as.double(K.hat),
              PACKAGE="simone")
    
    K.hat <- matrix(out$K.hat, ncol=p)
    # Re-normalize the solution to get thje associated precision matrix #
    K.hat <- Symmetrize(- diag(1/diag(S)) %*% (K.hat * upper.tri(K.hat))
                        + diag(1/diag(S)) )
  }
  
  ## SOLVING PENALIZED LIKELIHOOD WTH GLASSO ##
  else if ("glasso" %in% method) {

    ## Initializing additional matrices #
    if (is.null(Sigma.hat)) {
      Sigma.hat <- matrix(S + diag(diag(Rho)),nrow=p, ncol=p)
    }

    ## Test the positive definiteness of the starting matrix
    if (sum(eigen(Sigma.hat)$values < 0) == 0) {
      Beta  <- matrix(0, nrow=p-1, ncol=p)

      out <- .C("GLasso",
                as.integer (p),
                as.double  (S),
                as.double  (Rho),
                as.double  (eps),
                as.integer (maxIt),
                Sigma.hat = as.double(Sigma.hat),
                Beta 	= as.double(Beta),
                finalIt   = as.integer(1),
                PACKAGE   ="simone")
      
      Sigma.hat <- matrix(out$Sigma.hat, ncol=p)
      Beta <- out$Beta
      
      ## Computing K.hat by blockwise invertion of W ...
      out <- .C("Inverse",
                as.integer (p),
                as.double  (Sigma.hat),
                as.double  (Beta),
                K.hat   = as.double(K.hat),
                PACKAGE    = "simone")
      
      K.hat <- matrix(out$K.hat, ncol=p)
    } else {
      cat("\nWARNING! Penalty to small to start the algorithm with positive definite matrix - ABORTING AND RETURNING NULL\n")
      K.hat <- NULL
    }

  }

  ## Put names to the inferred Matrix K.hat if applicable ##
  if (!is.null(K.hat)) {
    dimnames(K.hat) <- list(colnames(data), colnames(data))
  }
  return(list(K.hat=K.hat, Sigma.hat=Sigma.hat))
}
