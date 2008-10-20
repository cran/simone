\name{SimDataAffiliation}
\alias{SimDataAffiliation}

\title{Simulate an affiliation graph and its associated Gaussian sample}

\description{Simulate a  gaussian sample from a gaussian graphical network with a latent  structure}

\usage{
SimDataAffiliation (p, n, proba.in, proba.out, alpha, proba.dust=0) 
}

\arguments{
  \item{p}{
    Scalar : number of nodes in the graph.
  }
  
  \item{n}{
    Scalar : number of Gaussian samples to simulate.
  }
  \item{proba.in}{
    Scalar : graph intra-class probability of edge existence.
  }

  \item{proba.out}{
    Scalar : graph  inter-class probability of edge existence.
  }

  \item{alpha}{ Vector : prior relative probability of the graph's
    classes.  Probabilities are internally normalized once added the
    dustbin probability, so as dust.proba remains unchanged and alpha
    sums to 1. Also defines the number of classes by its length.
  }

  \item{proba.dust}{Scalar : Prior probability of dustbin nodes. Default $0$.
  }

}

\value{
  Returns a List comprising :
  \item{cl.theo}{A vector of node classification }
  \item{K.theo}{A \eqn{p\times p}{\code{p x p}} precision matrix}
  \item{Sigma.theo}{A \eqn{p\times p}{\code{p x p}} covariance matrix}
  \item{data}{A \eqn{n \times p}{\code{n x p}} matrix that contains a
  $n$ samples of a Gaussian vector} 
}

\details{This function first draws node class belonging for the \code{p}
  nodes using a multinomial law parametered by \code{alpha} and
  \code{proba.dust}. The graph is then built as follows: edge existence
  or absence are drawn from a mixture of binomial distributions,
  parametered by \code{proba.in}, \code{proba.out}, \code{proba.dust}
  according to the node class belonging. The edges are then set to -1 or
  1 with probability 1/2 while ensuring symmetry. The obtained matrix is
  normalized so as to have a "heavy diagonal" of unitary values,
  ensuring inversability. This constitutes the graph's precision matrix
  \code{K.theo}. The theoretic covariance matrix \code{Sigma.theo} is
  obtained by inversing \code{K.theo}.

  A simulated data set \code{data} is then generated as \code{n} sample
  of a Gaussian vector with covariance matrix \code{Sigma.theo}.
}

\examples{
library(simone)

## The data set parameters
p <- 100
n <- 200
proba.in  <- 0.25
proba.out <- 0.025
alpha <- c(0.5,0.3,0.2)

# Simulate the data set
X  <- SimDataAffiliation (p,n,proba.in,proba.out,alpha,proba.dust=0.2)

## Plot the graph and summary the data set
g <- Gplot(X$K.theo, X$cl.theo, main="A modular network")
summary(X$data)
}

\author{C. Ambroise, J. Chiquet and A. Smith}
\keyword{datagen}
\keyword{graphs}
