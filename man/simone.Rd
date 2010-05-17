\name{simone}
\alias{simone}

\title{SIMoNe algorithm for network inference}

\description{%
  The \code{simone} function offers an interface to infer networks based
  on partial correlation coefficients in various contexts and methods
  (steady-state data, time-course data, multiple sample setup,
  clustering prior)
}

\usage{
simone(X,
       type       = "steady-state",
       clustering = FALSE,
       tasks      = factor(rep(1, nrow(X))),
       control    = setOptions())
}

\arguments{%
  \item{X}{%
    a \eqn{n\times p}{n x p} matrix of data, typically \eqn{n}{n}
    expression levels associated to the same \eqn{p}{p} genes. Can also
    be a \code{data.frame} with \eqn{n}{n} entries, each column
    corresponding to a variable (a gene). Specifying \code{colnames} to
    \code{X} may be convenient in view of results analysis, since it
    will be used to annotate the plots. Note that this is the only
    required argument.
  }
  
  \item{type}{%
    a character string indicating the data specification (either
    \code{"steady-state"} or \code{"time-course"} data). Default is
    \code{"steady-state"}.
  }
  
  \item{clustering}{%
    a logical indicating if the network inference should be perfomed by
    penalizing the edges according to a latent clustering discovered
    during the network structure recovery. Default is \code{FALSE}.
  }

  \item{tasks}{%
    A factor with \eqn{n}{n} entries indicating the task belonging
    for each observation in the multiple sample framework. Default is
    \code{factor(rep(1, nrow(X)))}, that is, all observations come from a
    unique homogeneous sample.
  }

  \item{control}{%
    A list that is used to specify low-level options for the
    algorithm, defined through the \code{\link{setOptions}} function.
  }
}

\details{Any inference method available ("neighborhood selection",
  "graphical-Lasso", "VAR(1) inference" and "multitask learning" - see
  \code{\link{simone-package}}) relies on an   optimization problem under
  the general form
  \if{latex}{\deqn{%
      \hat{\mathbf{\Theta}}_\lambda = \arg \max_{\mathbf{\Theta}} \mathcal{L}(\mathbf{\Theta};\mathbf{X}) - \lambda \ \textrm{pen}_{\ell_1}(\mathbf{\Theta},\mathbf{Z}),
  }}
  \if{html}{%
    \out{
<center>
    &Theta;<sub>hat</sub> (&lambda;) = argmax<sub>&Theta;</sub> L(&Theta;; data) - &lambda; * pen<sub>l1</sub>(&Theta;, Z),
</center>
  }
}
\if{text}{\deqn{Theta.hat(lambda) = argmax_Theta L(Theta) - lambda * pen-l1(Theta,Z),
}}
where \eqn{\mathcal{L}}{L} is the log-likelihood of the model (pseudo
  log-likelihood for "neighborhood selection") and \if{latex}{\eqn{\lambda}} \if{html}{\out{&lambda;}} \if{text}{lambda}
  is a penalty parameter which controls the sparsity level of the
  network. The \eqn{p\times p}{p x p} matrix
  \if{latex}{\eqn{\mathbf{\Theta}}}\if{html}{\out{&Theta;}}\if{text}{Theta} describes the parameters (basically, the
  edges) of the model, while \eqn{\mathbf{Z}}{Z} represents a latent
  clustering which is also estimated when the argument \code{clustering}
  is set to \code{TRUE}.

 The model and the penalty function
  \if{latex}{\eqn{\mathrm{pen}_{\ell_1}}}\if{html}{\out{pen<sub>l1</sub>}} \if{text}{pen-l1} differ according to the context
  (steady-state/time-course data, multitask learning and its associated
  coupling effect). For further details on the models, please check the
  papers listed in the reference section of
  \code{\link{simone-package}}.

  The criterion displayed during a SIMoNe run is the value of the
  penalized likelihood for the current values of the estimor  \if{latex}{\eqn{\hat{\mathbf{\Theta}}_\lambda}}\if{html}{\out{&Theta;<sub>hat</sub>(&lambda;)}} \if{text}{Theta.hat(lambda)} corresponding to
  a given value of the overall penalty level \if{latex}{\eqn{\lambda}}\if{html}{\out{&lambda;}}\if{text}{lambda}.
 
The following information criteria are also computed for any value of
\if{latex}{\eqn{\lambda}}\if{html}{\out{&lambda;}}\if{text}{lambda}
and part of the output of \code{simone}. The BIC (\emph{Bayesian
  Information Criterion})

 \if{latex}{\deqn{%
     \textrm{BIC}(\lambda) =
     \mathcal{L}(\hat{\mathbf{\Theta}}_\lambda;\mathbf{X}) -
     \textrm{df}(\hat{\mathbf{\Theta}}_\lambda) \frac{\log(n)}{2},
 }}
   
   \if{html}{\out{
<center>
    BIC(&lambda;) = L(&Theta;<sub>hat</sub>(&lambda;); data) - df(&Theta;<sub>hat</sub>(&lambda;)) log(n)/2,
</center>
}}
\if{text}{\deqn{BIC(lambda) = L(Theta.hat(lambda);X) - df(Theta.hat(lambda)) log(n)/2,}}
  and the AIC (\emph{Akaike Information Criterion})
\if{latex}{\deqn{%
    \textrm{AIC}(\lambda) = \mathcal{L}(\hat{\mathbf{\Theta}}_\lambda;\mathbf{X}) - \textrm{df}(\hat{\mathbf{\Theta}}_\lambda).
  }}  
  \if{html}{\out{
<center>
    AIC(&lambda;) = L(&Theta;<sub>hat</sub>(&lambda;); data) - df(&Theta;<sub>hat</sub>(&lambda;)) .
</center>
}
}
\if{text}{\deqn{AIC(lambda) = L(Theta.hat(lambda);X) - df(Theta.hat(lambda)).}}

}

\value{Returns an object of class \code{simone}, which is list-like and
  contains the following:

  \item{networks}{%
    a list with all the inferred networks stocked as adjacency matrices (the successive values of
    \if{latex}{\eqn{\hat{\mathbf{\Theta}}_\lambda}}\if{html}{\out{&Theta;}}\if{text}{Theta}
    controled by the penalty level
    \if{latex}{\eqn{\lambda}}\if{html}{\out{&lambda;}}\if{text}{lambda}). In
    the multiple sample setup, each element of the list is a list with
    as many entries as samples or levels in \code{tasks}.
  }

  \item{penalties}{%
    a vector of the same length as \code{networks}, containing the
    successive values of the penalty level.
  }

  \item{n.edges}{%
    a vector of the same length as \code{networks}, containing the
    successive numbers of edges in the inferred networks. In the multiple
    sample setup, \code{n.edges} is a matrix with as many columns as
    levels in \code{tasks}.
  }
  
  \item{BIC}{%
    a vector of the same length as \code{networks}, containing the value
    of the BIC for the successively estimated networks.    
  }

  \item{AIC}{%
    a vector of the same length as \code{networks}, containing the value
    of the AIC for the successively estimated networks.    
  }
 
  \item{clusters}{a size-\eqn{p}{p} factor indicating the class
    of each variable.
  }

  \item{weights}{%
    a \eqn{p\times p}{pxp} matrix of weigths used to adapt the penalty to
    each entry of the \code{Theta} matrix. It is inferred through the
    algorithm according to the latent clustering of the network.  When
    \code{clustering} is set to \code{FALSE}, all the weights are equal
    to "1", which mean no adaptive penalization.
  }

  \item{control}{%
    a list describing all the posterior values of the parameters used by
    the algorithm, to compare with the one set by the
    \code{\link{setOptions}} function. As a matter of fact, many of the options
    are defined depending on the nature of the data and can be automatically corrected during internal checks of the coherence of desired options to the characteristics of the data.
  }
  
}

\note{%
  If nothing particular is specified about the penalty through the
  \code{control} list (see \code{\link{setOptions}}), the default is to
  start from a value of
  \if{latex}{\eqn{\lambda}}\if{html}{\out{&lambda;}}\if{text}{lambda}
  that ensures an empty network. Then
  \if{latex}{\eqn{\lambda}}\if{html}{\out{&lambda;}}\if{text}{lambda} is
  progressively shrinked, as close to zero as possible. Along the
  shrinkage of
  \if{latex}{\eqn{\lambda}}\if{html}{\out{&lambda;}}\if{text}{lambda},
  only networks with different numbers of edges are kept in the final
  output.
}

\author{J. Chiquet}

\seealso{
  \code{\link{setOptions}}, \code{\link{plot.simone}},
  \code{\link{cancer}} and  \code{demo(package="simone")}. 
}
\examples{
## load the breast cancer data set
data(cancer)
attach(cancer)

## launch simone with the default parameters and plot results
plot(simone(expr))

## try with clustering now (clustering is achieved on a 30-edges network)
plot(simone(expr, clustering=TRUE, control=setOptions(clusters.crit=30)))

## try the multiple sample
plot(simone(expr, tasks=status))

detach(cancer)
}

\keyword{htest}
