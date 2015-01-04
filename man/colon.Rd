\name{colon}
\alias{colon}

\title{simplified gene expression data from Alon et al. (1999)}
\usage{
data(colon)
}

\description{
Gene expression data (2000 genes for 62 samples) from a DNA
 microarray experiments of colon tissue samples (Alon et al., 1999). 
}

\details{ This data set contains 62 colon tissue samples with 2000 gene expression levels. Among 62 samples, 40 are tumor tissues (coded 1) and 22 are normal tissues (coded -1).  
}

\value{
A list with the following elements:
  \item{x}{A matrix of 2000 columns and 62 rows standing for 2000 gene expression levels and 62 colon tissue samples. Each row corresponds to a patient.}
  \item{y}{A numeric vector of length 62 representing the tissue type (1 for tumor; -1 for normal).}
}

\source{The data were introduced in Alon et al. (1999).
\url{http://microarray.princeton.edu/oncology/affydata/index.html}}

\references{
Alon, U., Barkai, N., Notterman, D.A., Gish, K., Ybarra, S., Mack, D., and Levine, A.J. (1999).
``Broad patterns of gene expression revealed by clustering analysis of tumor and normal colon tissues probed by oligonucleotide arrays'',
\emph{Proceedings of the National Academy of Sciences}, \bold{96}(12), 6745--6750.\cr
}


\examples{
# load sdwd library
library(sdwd)

# load data set
data(colon)

# how many samples and how many predictors?
dim(colon$x)

# how many samples of class -1 and 1 respectively?
sum(colon$y == -1)
sum(colon$y == 1)
}


\keyword{datasets}
