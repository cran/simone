\name{Penalty}
\alias{Penalty}
\title{Penalty matrix constructor}

\description{Build a penalty matrix that can be adapted to a graph
  with affilation structure in case of available node classification. If
  no classification is available, a uniformly designed matrix is given.}

\usage{
  Penalty(X, ...)
}
\arguments{
  \item{X}{a \eqn{n\times p} matrix of data, taken from a
    i.i.d. size-$n$ sample.} 
  \item{...}{ Additional arguments are available, see Details}
}
\details{
  The \code{Penalty} function construct a penalty matrix according to
  the arguments passed by the user.

  Additional arguments are :
  \describe{
    \item{\code{rho}}{Base penalty value to use. Automatical yet
      conservative value is calculated by default}.
    \item{\code{risk}}{Probability of misclassification required
      for the calculation of \code{rho} when none is specified. Default
  \code{0.05}.}
    \item{classes}{Node classification for construction of a classified
      penalty matrix. Default \code{NULL}.}
    \item{multipliers}{\code{list} containing the 
      multiplicative coefficients for adjusting the base penalty \code{rho}
      applied to the edges, according to the involved nodes'
      classes.

      Default is \code{multipliers=list(intra=1,inter=1.05,dust=1.1)}}
  }
  
  If \code{rho} is not specified, a default value is computed, function
  of the data variance, the \code{risk} and the size \eqn{n} of the
  sample. This default value is very conservative and leads to very
  sparse graphs. As a matter of fact, the probability of
  misclassification holds for the $p^2$ potential edges. Thus, a
  correction can be applied by multiplying the risk by the expected
  number of edges. For instance, when one expects as many edges as
  nodes, a typical risk to apply is \eqn{0.05 /p^2 \times p = 0.05 \times p}.

  If \code{classes} is \code{NULL}, a uniform penalty matrix is
  returned. If a vector of classes belonging is specified, a classified
  version of the penalty matrix is built that enforces an affiliation
  structure, by penalizing more intra-class connections. The
  \code{multpliers} argument is a list with \code{inter}, \code{intra}
  and \code{dust} that permits to adjust the  penalty according to the
  involved nodes' classes.
}
\value{
  Return a \eqn{p\times p} matrix of penalties.
}

\author{C. Ambroise, J. Chiquet and A. Smith}

\seealso{\code{\link{SimDataAffiliation}},  \code{\link{Mplot}}}

\examples{
  library(simone)

  ## Data set generation
  p <- 100 # number of nodes
  n <- 200 # sample size
  proba.in  <- 0.15
  proba.out <- 0.005
  alpha <- c(.3,.2,.5)
  X <- SimDataAffiliation (p, n, proba.in, proba.out, alpha, proba.dust=0.2)

  ## Build a penalty matrix with an arbitrary base-value of rho and
  ## multiplier values that encourage affiliation structure
  P <- Penalty(X$data, rho = 0.1, classes = X$cl.theo,
               multipliers = list(intra=1,inter=1.5,dust=2))
	
  par(mfrow=c(1,2))
  Mplot(P, main="Penalty matrix")	
  Mplot(P, X$cl.theo, main="Ordered penalty matrix")
}

\keyword{cluster}
\keyword{graphs}

