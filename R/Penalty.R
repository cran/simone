## __________________________________________________________
##
## Penalty
## 
## Construct a Penalty matrix using various method, based
## upon an initial penalty value (scalar) rho
## various rules : Banerjee's, Meinshausen's and ours
##
## INPUT
##	- X : n*p data matrix
##
##   OPTIONAL INPUT
##      - rho: a scalar used as base penalty
##      - classes: a vector of nodes clusters
##      - multipliers: a vector of ...
##      - risk: probability of misclassification
##
## OUTPUT
##           - Rho: the penalty matrix
## __________________________________________________________

Penalty <- function(X, ...) {

  n <- nrow(X)
  p <- ncol(X)

  classes     <- sub.param("classes", NULL , ...)	
  risk        <- sub.param("risk",    0.05 , ...)	
  multipliers <- sub.param("multipliers", list(intra=1,inter=1.05,dust=1.1), ...)
  rho         <- sub.param("rho",     2/(n*BaseLambdaValue(X,risk=risk)), ...)	
  
  ## CASE 1: NO CLASS INFORMATION
  ##
  ## if no class information is available or if the nodes are spread
  ## over just one class, construct a uniform penalty
  ## ----------------------------------------------------------------
  if (nlevels(classes) <= 1) {
    Rho <- matrix(rho,p,p)

  ## CASE 2: USE CLASS INFORMATION
  ##
  ## if no class information is available, or if the nodes are spread
  ## over just one class construct a uniform penalty
  ## ----------------------------------------------------------------    
  } else {
    Rho <- matrix(rho * multipliers$inter, p, p)    
    for (k in levels(classes) ) {
      quadrant <- (classes == k)
      if (k=="D") {
        Rho[quadrant, ] <- rho * multipliers$dust
        Rho[, quadrant] <- rho * multipliers$dust
      } else {
        Rho[quadrant, quadrant] <- rho * multipliers$intra
      }
    }
  }
  
  ## Enforce the diagonal penalty to lambda0, to produce a
  ## invertible estimator of the covariance matrix
  if (n < p) {
    lambda0 <-  rep(max(colSums(var(X))),p)
    diag(Rho) <- lambda0 + lambda0/100
  }
  
  return(Rho)
} 

## __________________________________________________________
##
##	BaseLambdaValue
## 
## Calculates a penalty value according to the Data set X with
## the rule given in the paper (Matias' rule)
##
## INPUT
##	- X : n*p data matrix
##	- risk : probability of misclassification
##      - classes: vector of nodes clusters
##
## OUTPUT : if no class information is available or if nodes are
## spread over jsut one class, the output is a scalar lambda
## otherwise, a Q x Q matrix lambda_ql is sent back (Q is the number
## of classes)
## __________________________________________________________
##
BaseLambdaValue <- function(X,risk=0.05,classes=NULL) {

  lambda.dust <- -1
   
  n <- nrow(X)  
  p <- ncol(X)
  
  sigma <- sqrt(diag(var(X)))
  tnm2 <- -qt(risk/(2*p^2),df=n-2)

  ## CASE 1: NO CLASS INFORMATION
  ##
  ## If one class available or no class information available,
  ## compute the same lambda for everyone
  ## ----------------------------------------------------------------    
  if (nlevels(classes) <= 1) {
    S <- sigma %*% t(sigma)
    sig.i.sig.j <- S * upper.tri(S,diag=FALSE)
    lambda <- 2/n * sqrt(n-2 + tnm2^2) / (max(sig.i.sig.j) * tnm2  )

  ## CASE 2: USE CLASS INFORMATION
  ##
  ## if no class information is available, or if the nodes are spread
  ## over just one class construct a uniform penalty
  ## ----------------------------------------------------------------    
  } else {
    lambda <- matrix(lambda.dust,nlevels(classes),nlevels(classes))
    
    for (q in levels(classes)) {
      for (l in levels(classes)) {
        if ((q != "D") & (l != "D")) {
          S <- sigma %*% t(sigma)
          S[classes!=as.numeric(q),classes!=as.numeric(l)] <- 0
          sig.i.sig.j <- S * upper.tri(S,diag=FALSE)
          lambda[as.numeric(q),as.numeric(l)] <- 2/n * sqrt(n-2 + tnm2^2) / (max(sig.i.sig.j) * tnm2 )
        }
      }
    }
  }
  
  return(lambda)
}
