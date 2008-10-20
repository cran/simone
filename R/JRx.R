## __________________________________________________________
##
## JRx
##
## Authors : A.SMITH, J. CHIQUET (R)
##	     G.GRASSEAU (C)
##
## INPUT :
##	Tau   : node classification indicator matrix
##	X     : graph precision matrix
##	param : list of MixNetLaplace model parameters :
##		mu     : Q*Q matrix of class means
##		lambda : Q*Q matrix of class dispersions
##	 	alpha  : Q-vector of class proportions
##	C     : boolean : if true, use C code, else use R code
## OUTPUT :
##	scalar value
##
## Calculates the lower bound of the likelihood of interest
## for our MixNetLaplace variational algorithm
## __________________________________________________________
##
JRx <- function (Tau, X, param, C=TRUE) {

  if (is.null(param)) {
    cat("\n JRx : WARNING : null param argument, returning NULL !!!\n")
    return(NULL)
  }

  n      <- dim(Tau)[1]
  Q      <- dim(Tau)[2]

  mu 		<- param$mu
  lambda	<- param$lambda
  alpha		<- param$alpha

  u.tri  <- upper.tri(X)

  if (C==TRUE) {

    ret <- .C("JRx_C",
                                      	## Input
          as.integer (n),    		#  Number of Nodes
          as.integer (Q),     		#  Number of Classes
          as.double  (Tau),            	#  Tau(n,Q) : probas of node class belonging
          as.double  (X),             	#  X(n,n)   : adjacency matrix
          as.double  (mu),     		#  mu(Q,Q)  : mean matrix
          as.double  (lambda), 		#  lambda(Q,Q) : dispersion matrix
          as.double  (alpha),  		#  alpha(Q)    : class proportions
                                      	## Output
          J    = as.double(1),         	#  J return value  
          PACKAGE = "simone"
          )
    J <- ret$J

  } else {

    ## Mathématiquement, doit augmenter au cours des itérations

    ## J = - sum/iq(Tau.iq*ln(Tau.iq))
    ##     + sum/iq(Tau.iq*ln(alpha[q]))
    ##     + sum/i<j(sum/q,l(Tau.iq*Tau.jl*ln(fLaplace(...)))) 

    logTau <- pmax(log(Tau),log(.Machine$double.xmin))
    logalpha <- pmax(log(alpha),log(.Machine$double.xmin))
  
    ## Les 2 premiers termes de J...
    J <- - sum( Tau * logTau ) + sum ( Tau %*% logalpha )

    eps <- .Machine$double.xmin
    ## ... et le 3ème
    for (q in 1:Q) {
    
      for (l in 1:Q) {
      
          lnf.ijql <- lnFLaplace(X, mu[q,l], lambda[q,l]) 
          tau.ijql <- Tau[,q] %*% t(Tau[,l])
        
          J <- J + sum( tau.ijql[u.tri] * lnf.ijql[u.tri] )

      } # next l

    } # next q

  } # end if R version

  return(J)
}

