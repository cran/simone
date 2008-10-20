## __________________________________________________________
##
##                     IC Mstep
##
## INPUT
##	X : graph precision matrix
##	Tau : node classification indicator matrix
##	verbose : verbose mode
##
## OUTPUT
##   a list that contains
##	alpha  : Q-vector of class proportions
##	mu     : Q*Q matrix of class means
##	lambda : Q*Q matrix of class dispersions
##
## Internal
## M step for variational MixNetLaplace algorithm
## __________________________________________________________
##
IC.Mstep <- function(X, Tau, param.default=NULL, verbose=FALSE) {
  
  if (verbose==TRUE) {
    cat("\n     - Parameters estimation...")
  }

  ## INTIALIZE
  Q       <- ncol(Tau) # num classes
  n       <- nrow(Tau) # num nodes
  mu         <- matrix(0,Q,Q)
  lambda     <- matrix(0,Q,Q)
  theta.ijql <- matrix(0,n,n)
  u.trig     <- upper.tri(X)
  nodiag     <- upper.tri(X) | lower.tri(X)
  
  ## ESTIMATOR OF ALPHA
  if (is.null(param.default$alpha)) {
    alpha <- apply(Tau,2,mean)
    alpha <- pmax( alpha, .Machine$double.xmin )    
  } else {
    ## Supervized mode
    alpha <- param.default$alpha
  }
  
  for (q in 1:Q) {
    
    for (l in q:Q) {

      if (q != l) {
        ## calculation intermediary
        theta.ijql <- Tau[,q] %*% t(Tau[,l])
      
        ## ESTIMATOR OF MU is a weighted median
        
        if (is.null(param.default$mu)) {
          mu[q,l] <- weighted.median(X[nodiag],theta.ijql[nodiag])
        }
        
        ## ESTIMATOR OF LAMBDA is a weighted dispersion
        if (is.null(param.default$lambda)) {
          lambda[q,l] <- sum( theta.ijql[nodiag] * abs(X[nodiag] - mu[q,l] )) / sum(theta.ijql)
        }

      } else {

        ## calculation intermediary
        theta.ijql <- Tau[,q] %*% t(Tau[,l])
    
        ## ESTIMATOR OF MU is a weighted median
        if (is.null(param.default$mu)) {
          mu[q,l] <- weighted.median(X[u.trig],theta.ijql[u.trig])
        }
        
        ## ESTIMATOR OF LAMBDA is a weighted dispersion
        if (is.null(param.default$lambda)) {
          lambda[q,l] <- sum( theta.ijql[u.trig] * abs(X[u.trig] - mu[q,l]) ) / sum(theta.ijql)
        }
      }
        
    } # next l
    
  } # next q
  
  if (is.null(param.default$lambda)) {
    lambda[lower.tri(lambda)] <- lambda[upper.tri(lambda)]
  } else {
    lambda <- param.default$lambda
  }

  if (is.null(param.default$mu)) {
    mu[lower.tri(mu)] <- mu[upper.tri(mu)]
  } else {
    mu <- param.default$mu
  }
   
  return(list(mu = mu, lambda = lambda, alpha  = alpha))  
}

## __________________________________________________________
##
##                         Estep
##
## Authors : C. AMBROISE, J. CHIQUET, A. SMITH
##
## INPUT :
##	X : graph precision matrix
##	Tau : node classification indicator matrix
##	parameters : list :
##		alpha  : Q-vector of class proportions
##		mu     : Q*Q matrix of class means
##		lambda : Q*Q matrix of class dispersions
##	FP.maxIt : max number of iterations
##	eps   : convergence threshold
##	verbose : boolean : print to screen E step progress?
##
## OUTPUT : list :
##	Tau : node class probability matrix
##
## Internal
## E step for variational MixNetLaplace algorithm
## __________________________________________________________
##
IC.Estep <- function(	X, 
			Tau, 
			parameters, 
			FP.maxIt	= 50,
                  	eps		= 1e-6, 
			verbose		= FALSE ) {

  if (verbose==TRUE) {
    cat("\n     - Fixed point resolution... ")
  }
  
  ## =====================================
  ## 		INITIALIZING
  n <- nrow(Tau) # num nodes
  Q <- ncol(Tau) # num classes
  lambda <- parameters$lambda
  mu     <- parameters$mu
  alpha  <- parameters$alpha 

  ones     <- cbind(rep(1, n))
  class.names <- colnames(Tau)
 
  ## =====================================
  ##   SOLVING TAU WITH FIXED POINT ALGO

  ## convergence setup
  J 	  <- -Inf # JRx criterium
  E.Delta  <- Inf # convergence deltas
  E.Deltas <- c()
  E.It    <- 0 # iterations   
  convergence <- list( JRx.pb=FALSE, Iteration.pb=FALSE, Converged=FALSE, Pb=FALSE)

  ## convergence loop
  if (verbose) { cat(" iterate: ") }
  
  while ( (convergence$Pb==FALSE) && (convergence$Converged==FALSE) ) {

    E.It    <- E.It+1
    Tau.old <- Tau
    J.old   <- J
    
    if (verbose) { cat(" ",E.It) }
    
    ## prep current estimate of log(Tau)
    logTau <- matrix(0,n,Q)

    for (q in 1:Q) {
      for (l in 1:Q) { 
        
        ## normal laplacians
        Beta.ijql <- lnFLaplace(X, mu[q,l], lambda[q,l])
        diag(Beta.ijql) <- 0
        
        logTau[,q] <- logTau[,q] + apply(((ones %*% Tau[,l]) * Beta.ijql), 1, sum )

      } # next l
    } # next q
    
    ## Normalizing in the log space to avoid numerical problems
    logTau <- logTau - apply(logTau,1,max)

    ## Now going back to exponential with the same normalization
    Tau <- (matrix(1,n,1) %*% alpha) * exp(logTau)
    Tau <- pmin(Tau,.Machine$double.xmax)
    Tau <- pmax(Tau,.Machine$double.xmin)
    Tau <- Tau / (rowSums(Tau) %*% matrix(1,1,Q))
    Tau <- pmin(Tau,1-.Machine$double.xmin)
    Tau <- pmax(Tau,.Machine$double.xmin)   
    
    ## JRx criterium
    J <- JRx(Tau, X, parameters, C=TRUE)

     # convergence pb : JRx decreases
    convergence$JRx.pb <- (J < J.old)  	
    if (convergence$JRx.pb==TRUE) {
      Tau <- Tau.old
    }

    ## Delta criterium
    E.Delta <- J - J.old
    E.Deltas[E.It] <- E.Delta

    ## convergence ?
    convergence$Iteration.pb  <- (E.It    	> FP.maxIt) # have we hit iter.max ?
    convergence$Converged     <- (abs(E.Delta) 	< eps    ) # has the algo converged ?
    convergence$Pb <- (	convergence$Iteration.pb  || convergence$JRx.pb)

  } # repeat till convergence or itMax
 
  ## probabilistic class estimation
  colnames(Tau) <- class.names
  
  ## check non-convergence
  if ( convergence$Pb==TRUE ) {
    if (verbose) {  
      cat("   can't enhance the criteria anymore...\n" )
    }
  }

  return(Tau)
}	 

## __________________________________________________________
##
## 		InferClasses 
##
## Public
## InferVertexClasses algorithm
##
## INPUT :
##	K : graph precision matrix
##	Q : number of classes to estimate
##	note : classes may be lost
##
## OPTIONAL INPUT :
##	parameters	: list of model parameters for semi-supervised inference :
##		alpha, lambda, epsilon : see OUTPUT
##	init		: string : name of classification
##		initialisation tu use ("spectral","degree")
##
##	FP.maxIt	: max number of E step iterations
##	IC.maxIt	: max number of IC iterations
##	verbose		: boolean : print to screen estimation progress?
##
##	eps		: convergence thershold to use
##	degree.threshold : threshold under which node degrees are 
##		considered null ; null-degree nodes are considered
##		as dustbin nodes
##
## OUTPUT : list:
##	Tau    : node class indicator matrix
##	cl     : node classification vector (IC result)
##      mu     : Q*Q matrix of class means
##	lambda : Q*Q matrix of class dispersions
##	alpha  : Q-vector of class proportions
##	It     : final IC iteration attained
##	vJRx   : vector of JRx values (one per IC iteration)
##	Deltas : final convergence value
##	class.losses : boolean vector indicating at which iterations classes were lost
##
## CONVERGENCE FAILURE CODES :
##	A convergence failure code is a boolean vector giving information about
##	the algorithm's E step convergence. From left two right, the values denote :
##		JRx.pb : 	Calculations induced a decrease in the convergence criterium J(Rx)
##		Iteration.pb : 	No convergence before max. number of iterations reached
##		Converged :	True if the algorithm converged
##		Stabilized.pb :	The criterium has stabilized above the convergence threshold
##		Pb :		Is True if there is a problem.
##	Convergence failure codes are given at each E step when verbose=TRUE and when the convergence fails.
##
## Performs graph classification based on 
## a mixture of Laplace distributions
## __________________________________________________________
##
InferClasses <- function(K, Q, ...) {

  ## defaults for hidden optional arguments ##
  eps 			<- sub.param("eps"		, 1e-6	, ...)
  degree.threshold 	<- sub.param("degree.threshold"	, 1e-12	, ...)
  verbose 		<- sub.param("verbose"		, TRUE	, ...)
  param.default		<- sub.param("param.default"	, NULL	, ...)	
  FP.maxIt 		<- sub.param("FP.maxIt"		, 20	, ...)
  IC.maxIt 		<- sub.param("IC.maxIt"	        , 10	, ...)
  
  ## ==============================================
  ##            ALGORITHM INITIALISATION
  ## ==============================================
  
  if (verbose) { cat("   Clustering initialized with ")  }
  p <- ncol(K)
 
  if (Q<=1) {
    if (verbose) { cat(" Uniform ")  }
    # uniform initialization
    b.dust  <- (CalcDegrees(K)<degree.threshold) 
    init.cl <- rep("1", p)
    init.cl[b.dust] <- "D"
    K.d     <- K[!b.dust,!b.dust]
  } else { 
    if (verbose) { cat("Spectral clustering")  }
    ## default : spectral
    b.dust  <- (CalcDegrees(K)<degree.threshold) 
    K.d     <- K[!b.dust,!b.dust]
    init.cl <- rep("D", p)
    if (sum(b.dust) < p-Q) {
      init.cl[!b.dust] <- SpectralClustering(K.d,Q)
    }
  }

  ## If no connection...
  if (sum(b.dust) >= p-Q) {
    return(list(Tau 	     = matrix(1, 1, p), 
                cl 	     = rep("D", p),
                parameters   = NULL,
                iteration    = 0,
                J 	     = c(0),
                delta        = c(0),
                class.losses = c(0)) )
  }
  
  ## Removing the dust-bin
  Tau <- ClassNum2Indic(factor(init.cl[!b.dust]))
  cl  <- apply(Tau,1,which.max)
  
  ## ==============================================
  ##            EM-LIKE ALGORITHM
  ## ==============================================

  ## ______________________________________________
  ##
  ##		convergence setup
  ## ______________________________________________
  ##
  convergence <- list(JRx.pb	   = FALSE, 
                      iteration.pb = FALSE, 
                      Converged	   = FALSE, 
                      class.pb	   = FALSE,
                      Pb	   = FALSE)
  delta        <- c()
  class.losses <- c()
  iteration    <- 1

  ## ______________________________________________
  ##
  ##		1st M Step : 
  ##		OUTSIDE THE LOOP
  ## ______________________________________________
  ##
  parameters <- IC.Mstep(K.d, Tau, param.default)
  J <- JRx(Tau, K.d, parameters, C=TRUE)

  ## ______________________________________________
  ##
  ##		Tau convergence loop
  ## ______________________________________________
  ##
  while ((!convergence$Pb) && (!convergence$Converged)) {

    if (verbose) { cat("\n   IC iteration :",iteration)  }
    
    ## Keep backup of previous iteration
    Tau.old 	   <- Tau
    cl.old  	   <- cl
    parameters.old <- parameters
    
    ## -----------------------------------------------
    ##  		E Step
    ## -----------------------------------------------
    Tau <- IC.Estep(K.d, Tau.old, FP.maxIt = FP.maxIt,
                    parameters = parameters, verbose=verbose)    
    cl <- factor(apply(Tau, 1, which.max))
    
    ## Class loss check
    res <- CheckClassLoss(Tau, Q)
    class.losses <- c(class.losses, 1*res$loss)

    if (Q > res$Q) {
      ## don't allow class loss : return to previous estimations
      if (verbose) { cat("       Class loss at iteration", iteration) }

      Tau <- Tau.old
      cl  <- cl.old

      convergence$class.pb <- TRUE	# will stop algo

    } else {

      ## -----------------------------------------------
      ## 			M Step
      ## -----------------------------------------------
      parameters <- IC.Mstep(K.d, Tau, param.default, verbose=verbose)

      ## Calcul de la valeur du critère JRx
      J <- c(J, JRx(Tau, K.d, parameters, C=TRUE))
      
      ## Check JRx convergence
      if ( (J[iteration+1] < J[iteration]) && (iteration > 1) ){
        parameters  	   <- parameters.old
        convergence$JRx.pb <- TRUE
      }

    } # end if class.loss ... else ...

    ## -----------------------------------------------
    ## 			Convergence
    ## -----------------------------------------------
    delta <- c(delta, sum(abs(Tau.old-Tau)))

    convergence$Iteration.pb  <- (iteration > IC.maxIt) # have we hit iter.max?
    convergence$Converged     <- (delta[iteration] < eps) # has the algo converged?
    convergence$Pb <- (	convergence$Iteration.pb  || 
			convergence$JRx.pb  	  ||
			convergence$class.pb) # do we have a pb?

    # prepare next iter
    iteration <- iteration + 1

  } # repeat until convergence

  ## convergence warnings
  if (verbose==TRUE) {
    if(convergence$Pb==TRUE) {
      cat("   aborting\n")
    } else {
      cat("\n   IC converged in ", iteration-1," iterations\n") 
    }
  }

  # add dustbin class to classes from dusted estimation
  cl2 <- factor(rep("D", dim(K)[1]),levels=c(as.character(1:Q),"D"))

  if (sum(!b.dust)>0) {
    cl2[!b.dust] <- cl
    tmp <- matrix(0,p,Q)
    tmp[!b.dust,] <- Tau
    Tau2 <- data.frame(tmp,D=as.numeric((cl2=="D")))
    names(Tau2) <- levels(cl2)
    
    parameters$mu     <- cbind(parameters$mu,     rep(0, dim(parameters$mu)[1]))
    parameters$mu     <- rbind(parameters$mu,     rep(0, dim(parameters$mu)[2]))
    parameters$lambda <- cbind(parameters$lambda, rep(0, dim(parameters$lambda)[1]))
    parameters$lambda <- rbind(parameters$lambda, rep(0, dim(parameters$lambda)[2]))
  }
  
  return(list(Tau 	   = Tau2,
              cl 	   = factor(cl2),
              parameters   = parameters,
              iteration    = iteration,
              J 	   = J,
              delta        = delta,
              class.losses = class.losses))
} 

## lnFLaplace
##
## INPUT :
##	x : value(s)
##	mu : distribution mean
##	lambda : distribution dispersion
##
## Calculate the log of the laplacian distribution
## __________________________________________________________
##
lnFLaplace <- function (x, mu, lambda) {

  if (lambda > .Machine$double.xmin) {
    ## log-laplacian
    res <- -abs(x-mu)/lambda - log(2*lambda)
  } else {
    ## log-dirac : null where X_ij==mu_ql, -Inf elsewhere
    res <- 0 * (x == mu) + log(.Machine$double.xmin) * (x != mu)
  }
  
  res <- pmax(res,log(.Machine$double.xmin))
  res <- pmin(res,log(.Machine$double.xmax))
  
  return(res)
}

## __________________________________________________________
##
## weighted.median
##
## INPUT
##	x : vector of values
##	w : vector of weights
##
## OUTPUT
##	scalar value
##
## Calculate the weighted median of a vector
## __________________________________________________________
##
weighted.median <- function(x, w) { 
  
  ## null weights => uniform weight distribution 
  ##              => normal median
  s <- sum(w)
  if (s<.Machine$double.eps) {
    w <- rep(1, length(x)) 
    s <- sum(w)
  }

  ## The median is weighted by the distribution of w ## (the clusters)
  #x <- x[x!=0] 
  ind   <- order(x)
  kappa <- x[ind] 

  ## Fonction de répartition empirique de W = tauiq taujl
  fdrW <- cumsum(w[ind])/s
  
  # if there is one or less elements, beware...
  if (fdrW[1]>=0.5) {
    n <- 1
  } else {
    n <- max(which(fdrW<=0.5))
  }
  
  return(kappa[n])
} 

## __________________________________________________________
##
## CheckClassLoss
##
## INPUT
##	Tau : matrix of node class indicators
##	current.Q : current number of classes
##
## OUTPUT
##   a list that contains
##	ok : boolean vector : true where class proportions 
##		are greater or equal to 1/(number of nodes).
##		Where false, class is supposed lost.
##	loss : boolean : class loss?
##	Q : number of classes remaining
## __________________________________________________________
##
CheckClassLoss <- function(Tau, current.Q)  {

  alpha  <- colSums(Tau)
  ok     <- (alpha>=(1/dim(Tau)[1]))
  this.Q <- sum(ok)
  
  return( list(loss=(this.Q<current.Q), Q=this.Q, ok=ok ) )
}
