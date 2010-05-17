\name{simone-package}
\alias{simone-package}
\docType{package}
\title{Statistical Inference for MOdular NEtworks (SIMoNe)}

\description{%
  The \code{R} package \pkg{simone} implements the inference of
  co-expression networks based on partial correlation coefficients from
  either steady-state or time-course transcriptomic data. Note that with
  both type of data this package can deal with samples collected in
  different experimental conditions and therefore not identically
  distributed. In this particular case, multiple but related graphs are
  inferred at once.
  
  The underlying statistical tools enter the framework of Gaussian
  graphical models (GGM). Basically, the algorithm searches for a latent
  clustering of the network to drive the selection of edges through an
  adaptive \eqn{\ell_1}{l1}-penalization of the model likelihood.

  The available inference methods for edges selection and/or estimation
  include
  \describe{
    \item{neighborhood selection}{as in Meinshausen and Buhlman
      (2006), steady-state data only;
    }
    \item{graphical Lasso}{as in Banerjee \emph{et al}, 2008 and
      Friedman \emph{et al} (2008), steady-state data only;
    }
    \item{VAR(1) inference}{as in Charbonnier, Chiquet and
      Ambroise (2010), time-course data only;
    }
    \item{multitask learning}{as in Chiquet, Grandvalet and
      Ambroise (preprint), both time-course and steady-state data.
    } 
  }

  All the listed methods are \eqn{\ell_1}{l1}-norm based penalization,
  with an additional grouping effect for multitask learning (including
  three variants: "intertwined", "group-Lasso" and "cooperative-Lasso").

  The penalization of each individual edge may be weighted according to
  a latent clustering of the network, thus adapting the inference of the
  network to a particular topology. The clustering algorithm is
  performed by the \code{mixer} package, based upon Daudin, Picard and
  Robin (2008)'s Mixture Model for Random Graphs.

  Since the choice of the network sparsity level remains a current issue
  in the framework of sparse Gaussian network inference, the algorithm
  provides a full path of estimates starting from an empty network and
  adding edges as the penalty level progressively
  decreases. \emph{Bayesian Information Criteria} (BIC) and \emph{Akaike
  Information Criteria} (AIC) are adapted to the GGM context in order to
  help to choose one particular network among this path of solutions.

  Graphical tools are provided to summarize the results of a
  \code{simone} run and offer various representations for network
  plotting.
}

\details{
  Index:
  \preformatted{
cancer               Microarray data set for breast cancer
coNetwork            Random perturbations of a reference network
getNetwork           Network extraction from a SIMoNe run
plot.simone          Graphical representation of SIMoNe outputs
plot.simone.network  Graphical representation of a network 
rNetwork             Simulation of (clustered) Gaussian networks
rTranscriptData      Simulation of artificial transcriptomic data 
setOptions           Low-level options of the 'simone' function
simone               SIMoNe algorithm for network inference
  }
}

\section{Demos available}{%
  Beyond the examples of this manual, a good starting point is to have a
  look at the scripts available via \code{demo(package="simone")}. They
  make use of \code{simone}, main function in the package, in various
  contexts (steady-state or time-course data, multiple sample
  learning). All these scripts also illustrate the use of the different
  plot functions.

  \describe{

    \item{\code{demo(cancer_multitask)}}{%
      example on the \code{cancer} data set of the multitask approach
      with a cooperative-Lasso grouping effect across tasks. Patient
      responses to the chemiotherapy (pCR or not-pCR) split the data set
      into two distinct samples. Network inference is performed jointly
      on these samples and graphical comparison is made between the two
      networks.
    }
    
    \item{\code{demo(cancer_pooled)}}{%
      example on the \code{cancer} data set which is designed to compare
      network inference when a clustering prior is used or
      not. Graphical comparison between the two inferred networks
      (with/without clustering prior) illustrates how inference is
      driven to a particular network topology when clustering is
      relevant (here, an affiliation structure).
    }
    
    \item{\code{demo(check_glasso, echo=FALSE)}}{%
      example that basically checks the consistency between the
      \pkg{glasso} package of Friedman \emph{et al} and the \pkg{simone}
      package to solve the \eqn{\ell_1}{l1}-penalized Gaussian
      likelihood criterion suggested by Banerjee \emph{et al} in the
      \eqn{n>p}{n>p} settings. In the \eqn{n<p}{n<p} settings,
      \pkg{simone} provides sparser solutions than the \pkg{glasso}
      package since the underlying Lasso problems are solved with an
      active set algorithm instead of the shooting/pathwise coordinate
      algorithm.
    }
    
    \item{\code{demo(simone_multitask)}}{%
      example of multitask learning on simulated, steady-state data: two
      networks are generated by randomly perturbing a common ancestor
      with the \code{coNetwork} function. These two networks are then
      used to generate two multivariate Gaussian samples. Multitask
      learning is applied and a simple illustration of the use of the
      \code{setOptions} function is given.
    }

    \item{\code{demo(simone_steadyState)}}{%
      example of how to learn a single network from steady-state data. A
      sample is first generated with the \code{rNetwork} and
      \code{rTranscriptData} functions. Then the path of solutions of
      the \emph{neighborhood selection} method (default for single task
      steady-state data) is computed.
    }

    \item{\code{demo(simone_timeCourse)}}{%
      example of how to learn a single network from time-course data. A
      sample is first generated with the \code{rNetwork} and
      \code{rTranscriptData} functions and the path of solutions of the
      \emph{VAR(1) inference} method is computed, with and without
      clustering prior.
    }
  }
}

\author{\itemize{
    \item Julien Chiquet \email{julien.chiquet@genopole.cnrs.fr},
    \item Gilles Grasseau \email{gilles.grasseau@genopole.cnrs.fr},
    \item Camille Charbonnier \email{camille.charbonnier@genopole.cnrs.fr},
    \item Christophe Ambroise \email{christophe.ambroise@genopole.cnrs.fr}.
  }
}

\references{%
  J. Chiquet, Y. Grandvalet, and C. Ambroise (preprint). Inferring
   multiple graphical structures. \emph{preprint available on ArXiv}.
   \url{http://arxiv.org/abs/0912.4434}.

  C. Charbonnier, J. Chiquet,  and C. Ambroise (2010). Weighted-Lasso
  for Structured Network Inference from Time Course Data. \emph{Statistical
  Applications in Genetics and Molecular Biology}, vol. 9, iss. 1,
  article 15. \url{http://www.bepress.com/sagmb/vol9/iss1/art15/}
  
  C. Ambroise, J. Chiquet, and C. Matias (2009). Inferring sparse
  Gaussian graphical models with latent structure. \emph{Electronic Journal of
  Statistics}, vol. 3, pp. 205--238. \url{http://dx.doi.org/10.1214/08-EJS314}

  O. Banerjee, L. El Ghaoui, A. d'Aspremont (2008). Model Selection
  Through Sparse Maximum Likelihood Estimation. \emph{Journal of Machine
  Learning Research}, vol. 9, pp. 485--516.
  \url{http://www.jmlr.org/papers/volume9/banerjee08a/banerjee08a.pdf}

  J. Friedman, T. Hastie and R. Tibshirani (2008). Sparse inverse
  covariance estimation with the graphical Lasso. \emph{Biostatistics},
  vol. 9(3), pp. 432--441. \url{http://www-stat.stanford.edu/~tibs/ftp/graph.pdf}

  N. Meinshausen and P. Buhlmann (2006). High-dimensional graphs and variable
  selection with the Lasso. \emph{The Annals of Statistics}, vol. 34(3),
  pp. 1436--1462. \url{http://projecteuclid.org/DPubS/Repository/1.0/Disseminate?view=body&id=pdfview_1&handle=euclid.aos/1152540754} 

  J.-J. Daudin, F.Picard and S. Robin, S. (2008). Mixture model for
  random graphs. \emph{Statistics and Computing}, vol. 18(2),
  pp. 173--183. \url{http://www.springerlink.com/content/9v6846342mu82x42/fulltext.pdf}
  
}

\keyword{package}
