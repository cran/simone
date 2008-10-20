\name{Mplot}
\alias{Mplot}
\title{Graphics representation of a matrix}
\description{
  Display the image of a matrix which may be rearranged according to an optional
  node classification vector.
}

\usage{
  Mplot(A = NULL, cl = NULL, ...)
}

\arguments{
  \item{A}{The matrix to display. If \code{NULL}, only the color legend
  is displayed.}
  \item{cl}{Node classification vector. Default \code{NULL}.}
  \item{\dots}{See \emph{details} section.}
}

\details{Parameters can be specified in the call to \code{Mplot} to
  change the color palette, the axis labels, and general parameters.
  \itemize{
    \item \bold{Palette parameters}
    \describe{

      \item{\code{colors}}{Define the palette used when drawing the
	matrix. \code{rainbow}, \code{heat.colors},
	\code{terrain.colors}, \code{topo.colors}, \code{cm.colors} are
	available in \code{R}.  Two grey palettes (\code{grey} and
	\code{inv-grey}) and a light
	spectrum-like palette (\code{light}) are provided in this
	function. The number of color levels are chosen by adding ".x"
	where \code{x} is the level number. Default "light.64". }
      
     \item{\code{normalize}}{(default \code{TRUE}) \code{A} values are
	normalized so that generated colors belong to interval
	\code{1:length(colors)}. If \code{normalize = FALSE} then the
	normalization of the previous \code{Mplot} is used.
      }
	
      \item{\code{color.legend}}{(default \code{FALSE}) Display the
	color map.}
    }
     
    \item \bold{Label parameters}

    Values of the axis labels are set to the row and column matrix
    \code{dimnames} or indexes. Other parameters are:
  
   \describe{
     \item{\code{lab.col.cex}, \code{lab.row.cex}}{Expansion value of
       the column and rows labels. Default $1$.}
     \item{\code{lab.col.pos}, \code{lab.row.pos}}{Label orientation of the
       columns and rows: \code{"v"} for vertical or \code{"h"} for
       horizontal. Default "h".}
   } 

   \item \bold{Other parameters:} 
   \describe{
     \item{A.diag}{Boolean. If \code{TRUE}, diagonal terms are set to zero. }
     \item{cl.order}{Boolean. It true, the matrix is reordered to
       display classes in size ascending order.}
     \item{\code{margin}}{Margin between the matrix representation and
       the axis. Default $0$.}
     \item{\code{main}, \code{sub}}{Graphical options corresponding to
       the generic \code{plot} arguments. Default \code{NULL}.
     }
     \item{\dots}{Other arguments pass to functions \code{image} and  \code{axis}.
       See these functions for more.}
   }
 }
}
  
\value{
  No value is returned.
}
\author{G. Grasseau, A. Smith, C. Ambroise}

\seealso{\code{\link{Gplot}}, \code{\link{image}}}

\examples{
  library(simone)

  # Simulate a data set with obvious structure
  p <- 200
  n <- 1 
  proba.in  <- 0.2
  proba.out <- 0.01
  alpha <- c(.5,.3,.2)
  X <- SimDataAffiliation (p, n, proba.in, proba.out, alpha)

  # Create four pictures
  par(mfrow=c(2,2))
  Gplot(X$K.theo, X$cl.theo, main="Affiliation graph")
  Mplot(X$K.theo, main="Adjacency matrix", colors="light.64")
  Mplot(X$K.theo, X$cl.theo, lab.display=FALSE, main="Ordered matrix")
  Mplot(color.legend=TRUE)
}

\keyword{hplot}
\keyword{graphs}

