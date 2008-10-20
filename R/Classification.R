## __________________________________________________________
##
## ClassNum2Indic
##
## INPUT :
##	cl : vector of class labels
##
## OUTPUT : 
##	A matrix of corresponding class belonging indicators
##
## Transform class indexes to a matrix of class belonging indicators
## __________________________________________________________
##
ClassNum2Indic <- function (cl) {

  cl <- as.factor(cl)
  n  <- length(cl)

  x  <- matrix(0, n, length(levels(cl)))
  x[(1:n) + n * (unclass(cl) - 1)] <- 1

  dimnames(x) <- list(names(cl), levels(cl))

  return (x)
}

## __________________________________________________________
##
## SpectralClustering
##
## INPUT 
##	X : graph precision matrix
##	Q : number of classes to keep
##
## OUTPUT
##	node classification vector
##
## Spectral Clustering function (kmeans in eigen space)
## __________________________________________________________
##
SpectralClustering <- function (X,Q) {

  ## Taking absolute value and removing the diagonal
  W <- abs(X)
  diag(W) <- 0

  ## Normalized Laplacian
  D <- colSums(W)
  L <- diag(rep(1,ncol(W))) - diag(D^(-1/2)) %*% W %*% diag(D^(-1/2))
  
  ## go into eigenspace
  U <- eigen(L)$vectors

  if (Q > 1) {
    select <- (ncol(U)-Q+1):(ncol(U)-1)
    
    ## Applying the K-means in the eigenspace
    cl <- kmeans(U[,select], Q, nstart = 10, iter.max = 30)$cluster
  } else {
    cl <- as.factor(rep(1,ncol(W)))
  }
  
  return(cl)
}
