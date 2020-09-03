\name{sdwd-package}
\alias{sdwd-package}
\docType{package}
\title{
Sparse Distance Weighted Discrimination
}
\description{
This package implements the generalized coordinate descent (GCD) algorithm to efficiently compute the solution path of the sparse distance weighted discrimination (DWD) at a given fine grid of regularization parameters. Sparse distance weighted discrimination is a high-dimensional margin-based classifier.
}
\details{
\tabular{ll}{
Package: \tab sdwd\cr
Type: \tab Package\cr
Version: \tab 1.0.3\cr
Date: \tab 2020-02-16\cr
License: \tab GPL-2\cr
}

Suppose \code{x} is the predictors and \code{y} is the binary response. With a fixed value \code{lambda2}, the package produces the solution path of the sparse DWD over a grid of \code{lambda} values. The value of \code{lambda2} can be further tuned by cross-validation.

The package \code{\link{sdwd}} contains five main functions:\cr
\code{sdwd}\cr
\code{cv.sdwd}\cr
\code{coef.sdwd}\cr
\code{plot.sdwd}\cr
\code{plot.cv.sdwd}\cr
}
\author{Boxiang Wang and Hui Zou\cr
Maintainer: Boxiang Wang  \email{boxiang-wang@uiowa.edu}}
\references{
Wang, B. and Zou, H. (2016)
  ``Sparse Distance Weighted Discrimination", \emph{Journal of Computational and Graphical Statistics}, \bold{25}(3), 826--838.\cr
  \url{https://www.tandfonline.com/doi/full/10.1080/10618600.2015.1049700}\cr
  \cr
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software}, \bold{33}(1), 1--22.\cr
\url{https://www.jstatsoft.org/v33/i01/paper} \cr
\cr
Marron, J.S., Todd, M.J., Ahn, J. (2007)
``Distance-Weighted Discrimination"", 
  \emph{Journal of the American Statistical Association}, \bold{102}(408), 1267--1271.\cr
  \url{https://www.tandfonline.com/doi/abs/10.1198/016214507000001120}\cr
  \cr
Tibshirani, Robert., Bien, J., Friedman, J.,Hastie, T.,Simon,
  N.,Taylor, J., and Tibshirani, Ryan. (2012)
  Strong Rules for Discarding Predictors in Lasso-type Problems,
    \emph{Journal of the Royal Statistical Society, Series B}, \bold{74}(2), 245--266.\cr
  \url{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9868.2011.01004.x}\cr
  \cr
Yang, Y. and Zou, H. (2013)
``An Efficient Algorithm for Computing the HHSVM and Its Generalizations", 
  \emph{Journal of Computational and Graphical Statistics}, \bold{22}(2), 396--415.\cr
  \url{https://www.tandfonline.com/doi/full/10.1080/10618600.2012.680324}\cr
}
\keyword{DWD}
\keyword{SVM}
\keyword{high-dimensional}
\keyword{classification}
\keyword{package}

