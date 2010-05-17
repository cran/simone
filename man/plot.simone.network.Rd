\name{plot.simone.network}
\alias{plot.simone.network}

\title{Graphical representation of a network}

\description{%
  Displays the network contained in an object of class  \code{simone.network}.
}

\usage{%
\method{plot}{simone.network}(x,
     y    = NULL,
     type = "default", last.coord=FALSE, ...)}

\arguments{
  \item{x}{%
    an object of class \code{simone.network} to display. 
  }
  \item{y}{%
    an optional \code{simone.network} object to compare \code{x} with.
  }
  \item{type}{%
    network display types (see also details) are
    \describe{
      \item{\code{"circle"}}{displays the network nodes on a circle
	shape.}
      \item{\code{"circles"}}{displays the network nodes on circle
	shapes. The different circles correspond to node classes.
      }
      \item{\code{"cluster"}}{(default) displays network nodes (no
	underlying shape is used)}
      \item{\code{"overlap"}}{display a unique graph in which 2 graphs are overlaid}
      \item{\code{"4graphs"}}{(default) displays the two networks, the
	intersection and the symmetric difference between the two networks.}
    }
  }
    \item{last.coord}{%
    use last node coordinates if \code{TRUE}.
  }
  \item{...}{additionnal parameters}
}
\details{%
  This function plots a graph representation from a
  \code{simone.network} object. When available, the classification
  vector describing a partition of nodes is represented.
  
  Different node layouts (see \code{type} option) can be chosen to
  represent networks:
  \enumerate{%
    \item if a single \code{simone.network} object is provided, the
    available layouts are \code{cluster} (the default), \code{circle}
    (nodes are laid on one circle) and \code{circles} (nodes are laid on
    several circles, one circle for a node class);
    \item if two \code{simone.network} objects are provide, the available
    layouts are \code{4graphs} (the default, which displays both
    networks as well as the intersection and the difference between
    them) and \code{overlap} (which overlay two networks, representing common
    edges in gray, edges present in the first network in blue and edges
    present in the second network in red).
  }
}

\note{%
  When comparing two networks, the network with the more numerous edges
  should be passed as the first argument of \code{plot.simone.network}
  since the node positions for both networks will be computed so as the
  first graph is as readable as possible.
}

\author{G. Grasseau}

\seealso{\code{\link{plot.simone}}, \code{\link{simone}}}

\examples{
## data set and graph generation
lambda  <- 0.125
epsilon <- 0.00125
alpha <- c(1/3,1/3,1/3)

pi.affi <- matrix(epsilon,3,3)
diag(pi.affi) <- lambda

g1 <- rNetwork(p=200, pi=pi.affi, alpha=alpha)
g2 <- coNetwork(g1, delta=10)

plot(g1, type="cluster") # the default
plot(g1, type="circle" ) # one circle
plot(g1, type="circles" ) # one circle per cluster
plot(g1, g2, type="4graphs") # the default for multiple inputs
plot(g1, g2, type="overlap") # comparison of 2 networks on an unique graph
}

\keyword{hplot}
\keyword{graphs}
