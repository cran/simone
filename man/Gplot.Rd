\name{Gplot}
\alias{Gplot}

\title{Graph Plot}

\description{
  This function plots a  graph representation from its adjacency
  matrix using when available a   classification vector describing a
  partition of the nodes. .
}

\usage{
Gplot(A, cl = NULL, ...)
}

\arguments{
  \item{A}{
    graph adjacency matrix.
  }
  \item{cl}{
    node classification vector. Default \code{NULL}.
  }
  \item{...}{
    Additional optional arguments are available, see details.
  }
}

\value{
  A matrix of node coordinates is returned, along with a  graph plot.  
}

\details{This graph-drawing function is based on \code{C} routines
    extracted from the \code{sna} package and reorganized to adapt its
    possibilities in our framework (graphs are undirected, weighted and
    with no self-loops).

    Edges are coloured by sign (red for negative, blue for
    positive).  Edge widths are proportional to absolute edge weight,
    the maximal width being set to 5. Weak edges are dotted rather than
    solid.  The function returns the coordinates of the nodes which can
    be later passed back to Gplot.

    Additional arguments are :
    \describe{
    \item{\code{coord}}{
      A Matrix containing the 2D coordinates for the nodes (one node per row). Default \code{NULL}.
    }
    \item{\code{labels}}{
      A character vector containing the labels of the nodes. Default \code{NULL}.
    }
    \item{\code{main}, \code{sub}}{
      Graphical options corresponding to generic function \code{plot}
      arguments. Defaults \code{NULL}.
    }
    \item{\code{degree.threshold}}{
      Threshold under which edges are not plotted. Default \code{.Machine$double.xmin}.
    }
    \item{\code{display.isolate}}{
      Boolean to display the isolated nodes.  Default \code{TRUE}.
    }
    \item{\code{cols}}{
      Colour specification for each node cluster.
    }    
    \item{\code{max.edges}}{
      Maximal number of edges over
      which the graph will not be plot. Default \code{10000}.
    }
  }
}

\seealso{
  \code{\link{SimDataAffiliation}}, \code{\link{Mplot}}
}

\examples{
  library(simone)

  ## Graph generation
  p <- 100
  n <- 1
  proba.in  <- 0.15
  proba.out <- 0.005
  alpha <- c(.5,.3,.2)
  X <- SimDataAffiliation (p, n, proba.in, proba.out, alpha)

  ## Graphical representation
  g <- Gplot(X$K.theo, X$cl.theo, main = "A tri-modular graph")
}

\author{G. Grasseau, A. Smith and J. Chiquet}

\keyword{hplot}
\keyword{graphs}
