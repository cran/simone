\name{InferClasses}
\alias{InferClasses}

\title{ Vertex Classification Inference }

\description{
 Perform graph node classification using the graph's structure, based on 
 a latent structure model where edge weights  are Laplace distributed.
}

\usage{
  InferClasses(K, Q, ...)
}

\arguments{
  \item{K}{
    A \eqn{p \times p}{p x p} graph precision matrix.
  }

  \item{Q}{
    The number of classes to estimate.
  }
  \item{...}{ Additional arguments are available, see Details}
}

\value{
  Returns a list with the folowing components:
  \item{Tau}{
    A \eqn{p \times Q}{\code{p x Q}} node class probability matrix.
  }
  \item{cl}{
    Vector : node classification (IC result).
  }

  \item{parameters}{A list composed by the estimated parameters of the
    mixture of Laplace distributions: a \eqn{Q\times Q}{\code{QxQ}}
    matrix \code{lambda} of dispersions, a \eqn{Q\times Q}{\code{QxQ}}
    matrix \code{mu} of means and a size \eqn{Q}{\code{Q}} vector
    \code{alpha} of class proportions. If the function has been called
    with not-null \code{param.default} argument, \code{parameters} will
    be the same.  }
  
  \item{iteration}{
    Final iteration number
  }

  \item{J}{
    Vector of successive values of optimized criterion
  }

  \item{delta}{
    Vector of the successive convergence values.
  }

  \item{class.losses}{Boolean vector indicating if a class was lost for
  each iteration.}

}

\details{ \code{InferClasses} is based on the graph node classification
  algorithm MixNet, using a Laplacian mixture variant. The algorithm
  will stop if classes are lost during estimation.

  Additional arguments are :
  \describe{
    \item{\code{eps}}{
       Convergence threshold for algorithms. Default \code{1e-6}.
    }
    \item{\code{degree.threshold}}{
       Threshold under which edges are not considered. Default \code{1e-12}.
    }

    \item{\code{verbose}}{
      Boolean indicating if the progression is
      printed to screen. Default \code{TRUE}.
    }

    \item{\code{param.default}}{A \code{list} object for supervized
      learning, fixing the parameters in the mixture of Laplace
      distributions for edge weights: contains a \eqn{Q\times
      Q}{\code{QxQ}} matrix \code{lambda} of dispersions, a \eqn{Q\times
      Q}{\code{QxQ}} matrix \code{mu} of means  and a size
      \eqn{Q}{\code{Q}} vector \code{alpha} of class proportion. Default
      \code{NULL}.  }
    
    \item{\code{FP.maxIt}}{Maximum number of iteration
      for the fixed point resolution. Default \code{32}.
    }

    \item{\code{IC.maxIt}}{Maximum number of iteration for the full
      algorithm. Default \code{10}.
    }
  }
}

\references{
  J.-J. Daudin, F. Picard, et S. Robin. \emph{A mixture model for random
    graphs}. Stat. Comput., 18(2), p.173--183, 2008.
}

\seealso{
  \code{\link{SimDataAffiliation}}, \code{\link{Mplot}}, \code{\link{Gplot}}
}

\examples{
  library(simone)  

  ## Simulate a data set with obvious structure
  p <- 200
  n <- 1 
  proba.in  <- 0.2
  proba.out <- 0.01
  alpha <- c(.5,.3,.2)
  X <- SimDataAffiliation (p, n, proba.in, proba.out, alpha, proba.dust=0.2)

  ## estimate vertex classification
  res <- InferClasses(X$K.theo, 3)

  ## See results
  par(mfrow=c(2,2))
  Mplot(X$K.theo, X$cl.theo, main="Theoretical classification")
  Mplot(X$K.theo, res$cl, main="IC classification") 
  g1 <- Gplot(X$K.theo, X$cl.theo, main="Theoretical classification")
  g2 <- Gplot(X$K.theo, res$cl, main="IC classification", coord=g1) 
  print(table(X$cl.theo,res$cl))
}

\author{C. Ambroise, J. Chiquet and A. Smith }

\keyword{htest}
\keyword{cluster}
\keyword{graphs}
