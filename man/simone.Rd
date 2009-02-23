\name{simone}
\alias{simone}

\title{Statistical Inference methods for MOdular NEtworks (SIMoNe)}

\description{The \code{SIMoNe} algorithm iteratively combines edges
  estimation by using \eqn{\ell_1}-norm regularization method and node
  classification on the basis of a mixture model of Laplace
  distributions for edge weight. Intuitively, the node class information
  can be used to fine-tune the \eqn{\ell_1} penalty parameter for each edge,
  thus "focalising" more accurately on the correct graph structure.}

\usage{
res <- simone (X, Q, ...) 
}

\arguments{
  \item{X}{
     A \eqn{n \times p} empirical data matrix (e.g. gene epxressions).
  }

  \item{Q}{
    Scalar : max number of classes to estimate.
  }
  
  \item{...}{ Additional arguments are available, see Details.}
}

\value{Returns a \code{simone} object containing the results of the
  estimation. \code{simone} objects are list-like and contain the
  following :

  \item{cl}{Vector: the final estimation of node classification.}
  \item{K.hat}{Matrix:  the final estimation of the precision matrix
    (that is, the network adjacency matrix).}
  \item{K.hat.perfect}{contains the estimation of \code{InferEdges}
    using a classified penalty matrix built from the theoretic
    classification. Requires the additional argument \code{theo.cl}.
  }
  \item{delta}{A vector of convergence criterium values.}
  \item{iteration}{The final number of iterations.}
}

\details{By default, the algorithm alternates classes estimation and
  network estimation by himself, automatically constructing an penalty
  matrix that slighlty encourages  an affiliation strucure. The user can
  adapt the procedure by tuning the following parameters:
  \describe{
    \item{rho}{Base penalty value to use. Automatical yet
      conservative value is calculated by default.}

    \item{multipliers}{\code{list} containing the multiplicative
      coefficients for adjusting the base penalty \code{rho} applied to
      the edges, according to the involved nodes' classes.

      Default is \code{multipliers=list(intra=1,inter=1.05,dust=1.1)}
    }
    \item{\code{num.iter}}{
      The maximum number of iterations to attempt. Default \code{5}.
    }
    \item{\code{cl.theo}}{
      The theoretic node classification. Required for 'perfect'
     estimation. Default \code{NULL}.
    }
    \item{\code{verbose}}{
      Boolean to trace progress on screen. Default \code{TRUE}..
    }
  }
}

\references{
  C. Ambroise, J. Chiquet, and C. Matias. \emph{Penalized  maximum likelihood  inference for  sparse Gaussian
    graphical models with latent structure}. SSB report
  (20). INRA. Octobre 2008. \url{http://genome.jouy.inra.fr/ssb/preprint/}.
}

\seealso{
  \code{\link{InferEdges}},
  \code{\link{InferClasses}},
  \code{\link{Penalty}},
  \code{\link{SimDataAffiliation}},
  \code{\link{Gplot}}.
}

\examples{
library(simone)
## load the data set
data(cancer)

p <- length(gene.names)
cat("\n Breast cancer data set: real network with", p, "genes")

## SIMoNe parameters
Q <- 2
mult <- list(init=1,intra=1,inter=2,dust=4)
penalty.pcr <- Penalty(cancer.pcr, risk=0.1*p)[1,2]
penalty.not <- Penalty(cancer.notpcr, risk=0.1*p)[1,2]
rho.fracs <- c(1e-2,5e-2,7.5e-2,.25)
rhos.pcr  <- rho.fracs * penalty.pcr
rhos.not  <- rho.fracs * penalty.not

for (i in 1:length(rho.fracs)) {
  
  cat("\n Inferring network on the 'pCR' sample with rho =",rhos.pcr[i])
  res.pcr <- simone(cancer.pcr,Q,rho=rhos.pcr[i],multipliers=mult )
  cat("\n Inferring network on the 'not pCR' sample with rho =",rhos.not[i])
  res.not <- simone(cancer.notpcr,Q,rho=rhos.not[i],multipliers=mult)
  
  a <- readline(prompt="\nPress enter to continue...")
  
  layout(matrix(c(1,2),1,2),width=c(1,1))
  Gplot(res.pcr$K.hat, res.pcr$cl, cols=c(grey(0),grey(.8),grey(1)),
        display.labels=TRUE, labels=gene.names,main="pCR")
  Gplot(res.not$K.hat, res.not$cl, cols=c(grey(0),grey(.8),grey(1)),
        display.labels=TRUE,labels=gene.names,main="not pCR")
}
}

\author{ C. Ambroise, J. Chiquet and A. Smith }
\keyword{htest}
