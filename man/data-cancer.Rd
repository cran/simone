\name{cancer}
\alias{cancer}
\alias{cancer.full}
\alias{cancer.pcr}
\alias{cancer.notpcr}
\alias{gene.names}

\docType{data}

\title{Microarray data set for breast cancer} 

\usage{
  data(cancer)
}

\description{ This gene expression data set was provided by \emph{Hess
et al}. It concerns one hundred thirty-three patients with stage I--III
breast cancer.  Patients were treated with chemotherapy prior to
surgery. Patient response to the treatment can be classified as either a
pathologic complete response (pCR) or residual disease
(not-pCR). \emph{Hess et al} developed and tested a reliable multigene
predictor for treatment response on this data set, composed by a set of
26 genes having a high predictive value.

The data set is split into 2 parts (pCR and not pCR), on which network
inference algorithms should be applied independently: only individuals
from the same classes should be consider as independent and identically
distributed.}

\format{\code{cancer.full} is a data.frame of the complete data set (133
rows and 26 columns).
  
  \code{cancer.pcr} is a data.frame of pCR data (34 rows and 26
  columns).

  \code{cancer.notpcr} is a data.frame of not-pCR data (99 rows and 26
  columns).

  \code{gene.names} is a vector of names of the 26 selected genes.
}

\references{ 

Hess, K.R. and Anderson, K. and Symmans, W.F. and Valero, V. and
Ibrahim, N. and Mejia, J.A. and Booser, D. and Theriault, R.L. and
Buzdar, U. and Dempsey, P.J. and Rouzier, R. and Sneige, N. and Ross,
J.S. and Vidaurre, T. and Gomez, H.L. and Hortobagyi, G.N. and
Pustzai, L. \emph{Pharmacogenomic predictor of sensitivity to 
preoperative chemotherapy with Paclitaxel and Fluorouracil, Doxorubicin, and Cyclophosphamide in breast cancer}, Journal of Clinical Oncology, 2006, volume 24(26), pages 4236--4244.  

Natowicz, R. and Incitti, R. and Horta, E.G. and Charles, B. and
Guinot, P. and Yan, K. and Coutant, C. and Andre, F. and Pusztai,
L. Rouzier, R. \emph{Prediction of the outcome of a preoperative
chemotherapy in breast cancer using DNA probes that provide
information on both complete and incomplete response}, BMC
Bioinformatics, 2008, 9 (149).  
}


\examples{
## load simone library and the data set
library(simone)
data(cancer)

## summary of gene expression levels for pCR and not-pCR
summary(cancer.pcr)
summary(cancer.notpcr)
}

\keyword{datasets}
