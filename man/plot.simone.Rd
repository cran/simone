\name{plot.simone}
\alias{plot.simone}

\title{Graphical representation of SIMoNe outputs}

\description{Plots various outputs associated to a SIMoNe run.}

\usage{%
\method{plot}{simone}(x,
     output    = c("BIC", "AIC", "ROC", "PR", "path.edges",
                   "path.penalty", "sequence"),
     ref.graph = NULL,
     ask       = TRUE,  ...)
}

\arguments{
  \item{x}{
    output of a \code{simone} run (must be an object of class \code{simone})
  }

  \item{output}{%
    a vector of character string indicating which outputs must be plotted
    (picken from \code{"BIC"}, \code{"AIC"}, \code{"ROC"}, \code{"PR"},
    \code{"path.edges"}, \code{"path.penalty"} or
    \code{"sequence"}). Default is to plot everything possible.
  }
  
  \item{ref.graph}{%
    a network of reference provided through an adjacency matrix that is
    used to compute the ROC and PR curves. Only required if \code{"ROC"}
    and \code{"PR"} belongs to the \code{output} argument.
  }

  \item{ask}{%
    a logical indicating if the graphics device should be
    interactive. Default is \code{TRUE}.
  }

  \item{\dots}{
    Additonal arguments for generic \code{plot} (such as \code{main = "my
    title"}).
  }
}

\author{J. Chiquet}

\details{Here are some details about the plots possibly produced:
  \itemize{
    
    \item{%
      If \code{"BIC"} belongs to the \code{output} argument, a plot
      representing the Bayesian Information Criterion as a function of
      each network inferred by \code{simone} is displayed.
    }

    \item{%
      If \code{"AIC"} belongs to the \code{output} argument, a ploy
      representing the Akaike Information Criterion as a function of
      each network inferred by \code{simone} is displayed.
    }

    \item{%
      If \code{"ROC"} belongs to the \code{output} argument and
      \code{ref.graph} is specified, the ROC curve (Receiver
      Operating Characteristic) is plotted by representing true positive
      rate vs. false positive rate.
    }

    \item{%
      If \code{"PR"} belongs to the \code{output} argument and
      \code{ref.graph} is provided by the user, the PR curve
      (Precision/Recall) is plotted by representing positive predicted
      values vs. true positive rate.
    }

    \item{%
      If \code{"path.penalty"} belongs to the \code{output} argument, a
      regularization path is plotted by representing the value of each
      entry of the \code{Theta} matrix (that is, of each edge) vs. the
      penalty level \eqn{\lambda}{lambda}: there are as many values for
      the penalty as networks stocked in the \code{simone} object \code{x}.
    }

    \item{%
      If \code{"path.edges"} belongs to the \code{output} argument, a
      regularization path is plotted by representing the value of each
      entry of the \code{Theta} matrix (that is, of each edge) vs. the
      degree of freedom in \code{Theta} (that is, the number of edges in the
      current network). This is done for all the network stocked in the
      \code{simone} object \code{x}.
    }

    \item{%
      If \code{"sequence"} belongs to the \code{output} argument, an
      interactive plot is provided by starting from the empty network and
      adding the edges by successively covering the networks stocked in
      the \code{simone} object \code{x}.
    }    
  }  
}

\note{%
  If the user asked for \code{"PR"} and \code{"ROC"} curves yet did not specify a
  network of reference, these curves will not be plotted (no warning or
  error will be specified).
}

\seealso{\code{\link{simone}}.}

\examples{
## data set and network generation
g    <- rNetwork(p=50, pi=50)
data <- rTranscriptData(300,g)
attach(data)

## running simone
res <- simone(X, type="steady-state")

## plotting the results: just the ROC curve
plot(res, output=c("ROC"), ref.graph=g$A)

## plotting the results: just the path as a function of the penalty
plot(res, output=c("path.penalty"), main="I want to put my own title")

## plotting the results: everything possible (the default)
plot(res)

detach(data)
}

\keyword{hplot}
