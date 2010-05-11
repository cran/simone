\name{setOptions}
\alias{setOptions}

\title{Low-level options of a SIMoNe run}

\description{%
  This function is intended to design low-level uses of SIMoNe by
  specifying various parameters of the underlying algorithms.
}

\usage{%
setOptions(normalize      = TRUE,
           verbose        = TRUE,
           penalties      = NULL,
           penalty.min    = NULL,
           penalty.max    = NULL,
           n.penalties    = 100, 
           edges.max      = Inf,
           edges.sym.rule = NULL,
           edges.steady   = "neighborhood.selection",
           edges.coupling = "coopLasso",
           clusters.crit  = "BIC",
           clusters.meth  = "bayesian",
           clusters.qmin  = 2,
           clusters.qmax  = 4)
}

\arguments{
  \item{normalize}{%
    logical specifying wether the data should be normalized to unit
    variance. The normalization is made task-wisely in the multiple
    sample setting. Default is \code{TRUE}.
  }

  \item{verbose}{%
    a logical that indicates verbose mode to display
    progression. Default is \code{TRUE}.
  }
 
  \item{penalties}{%
    vector of \emph{decreasing} penalty levels for the network
    estimation. If \code{NULL} (the default), an appropriate vector will
    be generated in \code{simone} with \code{n.penalties} entries,
    starting from \code{penalty.max} and shrinked to
    \code{penalty.min}.
  }
  
  \item{penalty.min}{%
    The minimal value of the penalty that will be tried for network
    inference. If \code{NULL} (the default), it will be set in \code{simone} to
    \code{1e-5} for the monotask framework and to \code{1e-2} for the
    multitask framework.
  }

  \item{penalty.max}{%
    The maximal value of the penalty that will be tried for network
    inference. If \code{NULL} (the default), it will be set to a value
    that provoques an empty granph. Default is \code{NULL}.
  }

  \item{n.penalties}{%
    integer that indicates the number of penalties to put in the
    \code{penalties} vector. Default is \code{100}.
  }

  \item{edges.max}{%
    integer giving an upper bound for the number of edges to select: if
    a network is inferred along the algorithm with a number of edges
    overstepping \code{edges.max}, it will stop there. Default
    is \code{Inf}.
  }

  \item{edges.steady}{%
    a character string indicating the method to use for the network
    inference associated to steady-state data, one task
    framework. Either \code{"graphical.lasso"} or
    \code{"neighborhood.selection"}. Default is  the later.
  }

  \item{edges.coupling}{%
    character string (either \code{"coopLasso"}, \code{"groupLasso"}
    or \code{"intertwined"}) that indicates the coupling method across
    task in the multiple sample setup. Defautl is \code{"coopLasso"}.
  }
  
  \item{edges.sym.rule}{%
    character string (\code{"AND"}, \code{"OR"}, \code{"NO"}) for
    post-symmetrization of the infered networks. Enforced to \code{"NO"}
    for time-course data (directed network) and set to \code{"AND"} as
    default for steady-state data (undirected network).
  }

  \item{clusters.crit}{%
    criterion to select the network that is used to find an underlying
    clustering. Either \code{"BIC"}, \code{"AIC"} or an integer for the
    number of edges. Default is \code{"BIC"}.  
  }

  \item{clusters.qmin}{%
    minimum number of classes for clustering. Default is 2.
  }

  \item{clusters.qmax}{%
    maximum number of classes for clustering. Default is 4.
  }

  \item{clusters.meth}{%  
    character string indicating the strategy used for the estimation:
    \code{"variational"}, \code{"classification"}, or
    \code{"bayesian"}. See the \pkg{mixer} package for further
    details. Default is \code{"bayesian"}.
  }
}

\note{%
  If the user specifies its own \code{penalties} vector, all the
  networks inferred during the algorithm will be kept, even if they
  share the very same number of edges.

  On the other hand, if you only specify \code{penalty.max} and/or
  \code{penalty.min} and/or \code{n.penalties}, the algorithm will only
  kept the networks who show different numbers of edges. That is to say,
  the number of networks stocked in the output of \code{simone}
  generally does not have a length equal to \code{n.penalties}.
  
}

\value{
  A list that contains all the specified parameters.
}

\author{J. Chiquet}

\seealso{
  \code{\link{simone}}.
}
\examples{
## generate an object (list) with the default parameters
setOptions()
}

\keyword{misc}
