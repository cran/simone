\name{rNetwork}
\alias{rNetwork}

\title{Simulation of (clustered) Gaussian networks}

\description{Simulates a network with various structures.}
\usage{%
rNetwork(p,
         pi,
         alpha    = c(1),
         directed = FALSE, 
         name     = "a network",
         signed   = TRUE)
}

\arguments{%
  \item{p}{%
    the number of nodes of the simulated network.
  }
  \item{pi}{%
    a matrix of cluster connectivity (see details).
  }
  \item{alpha}{%
    a vector of cluster proportions.
  }
  \item{directed}{%
    a logical indicating the directedness of the network.
  }
  \item{name}{%
    a character string  indicating the name of the network.
  }
  \item{signed}{%
   a logical indicating whether partial correlations should be signed or
  all kept positive.
  }
}
\details{%
  Matrix \code{pi} should be a square matrix of the same size as vector
  \code{alpha}. When the network is not directed, \code{pi} should be
  symmetric. When the graph is directed, entry
  \ifelse{latex}{\eqn{\pi_{q\ell}}}{\out{&pi;}_ql}    corresponds to
  edges heading from class \eqn{q} to class  \eqn{\ell}{l}.

  Entries of \code{pi} can be either integers or real numbers. If they are
  integers, they are considered as the exact number of edges required from
  one class to another. Otherwise, they are considered as connectivity
  probabilities between classes. They should therefore sum up to at most
  1. If they do not sum up to one excatly, the remaining value is
  considered as the probability for a node to belong to the dust class
  (connected to no other node).
}

\value{%
  Returns an object of class \code{simone.network}, that is, a list
  comprising
  \item{A}{the \eqn{p \times p}{p x p} adjacency matrix of the
    network, filled with 0 and 1's, which is symmetric if
    \code{directed} is \code{FALSE}.}
  \item{Theta}{a \eqn{p \times p}{p x p} matrix of parameters of the
    associated Gaussian model, which depends on the directedness of the
    network: if directed, \code{Theta} contains the parameters of a
    VAR(1) model; if undirected, \code{Theta} is the concentration
    matrix (inverse of the covariance matrix) of a Gaussian vector}
  \item{directed}{a logicial indicating the directedness of the
    network.}
  \item{clusters}{a size-\eqn{p}{p} factor indicating the node
    class. The number of levels is determined by the number of columns
    of the matrix of connectivity \code{pi}: the levels are labeled
    \eqn{1,\dots,Q}{1,...,Q} where \eqn{Q}{Q} is the number of clusters.
  }
  \item{name}{a character string containing the name of the network.}
}

\author{J. Chiquet, C. Charbonnier}

\seealso{
  \code{\link{coNetwork}}, \code{\link{plot.simone.network}}.
}
\examples{
## generate an Erdos-Renyi network with 50 nodes and Pr of edges = 0.1
plot(rNetwork(p = 50, pi = 0.1, name = "an Erdos-Renyi network"))

## generate an network with 15 nodes and 25 randomly selected edges
plot(rNetwork(p = 15, pi = 25, name = "a 25 edges network"))

## generate an undirected network with an affiliation structure
PI <- matrix(c(15,2,2,50),2,2)
alpha <- c(1/3,2/3)
plot(rNetwork(p = 20, pi = PI, alpha = alpha,
                      name = "Affiliation, fixed num of edges"))

## generate a directed network with hubs
PI <- t(matrix(c(0.2,0.1,0.4,0,0.05,0.15,0,0.4,rep(0,8)),4,4))
alpha <- c(1/20,1/20,9/20,9/20)
plot(rNetwork(p = 55, pi = PI, alpha = alpha, directed = TRUE,
                      name = "Hubs structured network"))
}

\keyword{graphs}
\keyword{datagen}
