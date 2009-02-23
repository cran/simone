\name{InferEdges}
\alias{InferEdges}

\title{Edge Inference}

\description{Estimate the inverse covariance matrix from a
  i.i.d. size--\eqn{n}{n} sample of a multivariate normally distributed
  random vector.}

\usage{
  InferEdges(data, penalty, method="glasso", ...)
}


\arguments{ \item{data}{A \eqn{n \times p}{n x p} data matrix containing
    i.i.d. size--\eqn{n}{n} sample taken from a multivariate normally
    distributed random size--\eqn{p}{p} vector.
  }

  \item{penalty}{Penalty to use. Can be a numerical matrix of size
    \eqn{p \times p}{p x p} or a scalar value. If \code{NULL}, a default
    conservative penalty is calculated that will lead to very sparse
    graph.}

  \item{method}{ A string that defines the method to use for the
    estimation of the inverse covariance matrix: either \code{"glasso"},
    \code{"regressionAND"} or \code{"regressionOR"}.  Default is
    \code{"glasso"}.
  }

  \item{...}{ Additional arguments are available, see Details}
}

\value{
  Return a list with the two following components:

  \item{Sigma.hat}{
    The \eqn{p \times p}{\code{p x p}} estimated covariance matrix.
  }

 \item{K.hat}{ The \eqn{p \times p}{\code{p x p}} estimated inverse
   covariance (or precision) matrix.
 }

 Note that \code{Sigma.hat} is \code{NULL} for "regressionOR" and
 "regressionAND", since the precision matrix \code{K.hat} is the only
 one estimated with these methods.  }

\details{ 

  \code{InferEdges} is a wrapper for accessing our implementation in
  \code{C} of several algorithms for estimating inverse covariance
  matrices. Given the link between such matrices and corresponding graph
  precision matrices, we currently substitute the latter for the former,
  hence the name \code{InferEdges}. The implemented inference algorithms
  are :
  \describe{
    \item{The "glasso" method,}{which solves a
      \eqn{\ell_1}{l1}--penalized likelihood problem (Banerjee et al,
      2008) based upon the GLasso approach (see Friedman et al, 2007).}

    \item{The "regressionAND" method,}{which solves \eqn{p}{p}
      independent \eqn{\ell_1}{l1}--penalized regressions with an AND rule
      for symmetrization (see Meinshausen and B\"ulhman).}
    
    \item{The "Regression OR",}{which solves \eqn{p}{p} independent
      \eqn{\ell_1}{l1}--penalized regressions with an OR rule for
      symmetrization (see Meinshausen and B\"ulhman).}
  }
  
  The penalty term can be a scalar or a matrix. For the latter, the
  penalty is applied term-to-term to the inverse covariance matrix
  estimator, thus penalizing each entry differently.

  Additional arguments are :

  \describe{

    \item{\code{Sigma.hat}}{ \eqn{p \times p}{p x p} matrix. Starting point of
      the algorithm. If NULL, use \eqn{S_n+ diag(penalty)}{S_n+
	diag(penalty)} where \eqn{S_n}{S_n} is the empirical covariance
      matrix. Default \code{NULL}.
    }

    \item{\code{eps}}{ Scalar. Convergence threshold for the
      algorithms. Default \code{1e-12}.
    }
    
    \item{\code{maxIt}}{
       Maximum number of iterations for block-wise coordinate
       algorithm. Default \code{1e4}.
    }
  }
}

\references{
  Banerjee, O., El Ghaoui, L. and d'Aspremont, A. \emph{Model
  selection through sparse maximum likelihood estimation for
  multivariate Gaussian or binary data},  Jour. Mach. Learn. Res., 9,
  p.~485--516, 2008.

  Friedman, J., Hastie, T. and Tibshirani, R.  \emph{Sparse inverse
  covariance estimation with the graphical lasso}, 9(3), p.~ 432--441,
  Biostatistics, 2008.

  Meinshausen, N. and B\"uhlmann, P.  \emph{High-dimensional graphs and
    variable selection with the lasso}, Ann. Statist., 34(3),
  p.~1436--1462, 2006
}

\seealso{
  \code{\link{SimDataAffiliation}}, \code{\link{Gplot}}, \code{\link{Mplot}}
}

\examples{
  library(simone)

  ## Generating a graph with an associated Gaussian sample
  p <- 100
  n <- 200
  proba.in  <- 0.15
  proba.out <- 0.005
  alpha <- c(.6,.4)
  X <- SimDataAffiliation (p, n, proba.in, proba.out, alpha)

  ## Network inference
  rho <- 0.18
  res <- InferEdges(X$data, rho)

  ## Results, plotting and comparison
  par(mfrow=c(2,2))
  g <- Gplot(X$K.theo, X$cl.theo, main="Theoretical graph")
  Mplot(X$K.theo, X$cl.theo, main="Theoretical Mplot")
  Gplot(res$K.hat, coord=g, main="GLasso Inference")
  Mplot(res$K.hat, X$cl.theo, main="Inferred Mplot")
}

\author{J. Chiquet, based upon ealier work of J. Friedman and R. Tibshirani}

\keyword{htest}
\keyword{graphs}
