\name{sdwd-package}
\alias{sdwd-package}
\docType{package}
\title{
Sparse Distance Weighted Discrimination
}
\description{
  This package implements the generalized coordinate descent (GCD) algorithm to efficiently compute the solution path of the sparse Distance Weighted Discrimination (DWD) at a given fine grid of regularization parameters.
}
\details{
\tabular{ll}{
Package: \tab sdwd\cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2015-08-05\cr
License: \tab GPL-2\cr
}

Suppose \code{x} is the predictors and \code{y} is the binary response. With a fixed value \code{lambda2}, the package produces the solution path over a grid of \code{lambda} values.

The package \code{\link{sdwd}} contains five main functions:\cr
\code{sdwd}\cr
\code{coef.sdwd}\cr
\code{predict.sdwd}\cr
\code{print.sdwd}\cr
\code{plot.sdwd}\cr
}
\author{Boxiang Wang and Hui Zou\cr
Maintainer: Boxiang Wang  \email{boxiang@umn.edu}}
\references{
Wang, B. and Zou, H. (2015)
  ``Sparse Distance Weighted Discrimination", \emph{Journal of Computational and Graphical Statistics}, forthcoming. \cr
  \url{http://arxiv.org/abs/1501.06066}\cr
  \cr
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software}, \bold{33}(1), 1--22\cr
\url{http://www.jstatsoft.org/v33/i01/paper} \cr
\cr
Marron, J.S., Todd, M.J., Ahn, J. (2007)
``Distance-Weighted Discrimination"", 
  \emph{Journal of the American Statistical Association}, \bold{102}(408), 1267--1271\cr
  \url{https://faculty.franklin.uga.edu/jyahn/sites/faculty.franklin.uga.edu.jyahn/files/DWD3.pdf}\cr
  \cr
Tibshirani, Robert., Bien, J., Friedman, J.,Hastie, T.,Simon,
  N.,Taylor, J., and Tibshirani, Ryan. (2012)
  Strong Rules for Discarding Predictors in Lasso-type Problems,
    \emph{Journal of the Royal Statistical Society, Series B}, \bold{74}(2), 245--266\cr
  \url{http://statweb.stanford.edu/~tibs/ftp/strong.pdf}\cr
  \cr
Yang, Y. and Zou, H. (2013)
``An Efficient Algorithm for Computing the HHSVM and Its Generalizations", 
  \emph{Journal of Computational and Graphical Statistics}, \bold{22}(2), 396--415\cr
  \url{http://users.stat.umn.edu/~yiyang/resources/papers/JCGS_gcdnet.pdf}\cr
}
 
\keyword{DWD}
\keyword{SVM}
\keyword{high-dimensional}
\keyword{classification}
\keyword{package}

