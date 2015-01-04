\name{sdwd-internal}
\title{internal sdwd functions}
\alias{coef.sdwdNET}
\alias{cv.sdwdNET}
\alias{cvcompute}
\alias{err}
\alias{error.bars}
\alias{getmin}
\alias{getoutput}
\alias{lambda.interp}
\alias{lamfix}
\alias{nonzero}
\alias{sdwdpath}
\alias{zeromat}
\description{Internal sdwd functions.}
\usage{
cv.sdwdNET(outlist, lambda, x, y, foldid)
cvcompute(mat, foldid, nlams)
err(n, maxit, pmax)
error.bars(x, upper, lower, width=0.02, ...)
getmin(lambda, cvm, cvsd)
getoutput(fit, maxit, pmax, nvars, vnames)
lambda.interp(lambda, s)
lamfix(lam)
nonzero(beta, bystep=FALSE)
zeromat(nvars, nalam, vnames, stepnames)
}
\author{Boxiang Wang and Hui Zou\cr
Maintainer: Boxiang Wang  \email{boxiang@umn.edu}}
\references{Wang, B. and Zou, H. (2015)
  ``Sparse Distance Weighted Discrimination".%, \emph{University of Minnesoa Statistics Technical Report.}
  \cr
  \url{}\cr
Yang, Y. and Zou, H. (2013)
``An Efficient Algorithm for Computing the HHSVM and Its Generalizations", 
  \emph{Journal of Computational and Graphical Statistics}, \bold{22}(2), 396--415\cr
  \url{http://users.stat.umn.edu/~yiyang/resources/papers/JCGS_gcdnet.pdf}\cr
  \cr
  Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software}, \bold{33}(1), 1--22\cr
\url{http://www.jstatsoft.org/v33/i01/paper}
\cr
 }

\details{These are not intended for use by users. \code{coef.sdwdNET} computes the coefficient of the \code{sdwd} object. \code{cv.sdwdNET} does cross-validation for the \code{sdwd} object. \code{cvcompute} computes the mean and the standard deviation of the cross-validation error. \code{err} obtains the error message from fortran code. \code{error.bars} helps to plot the cross-validation error curve. \code{getmin} addresses the best lambda through the cross-validation either using or not using the one-standard-deviation rule. \code{getoutput} organizes the output of the \code{sdwd} object. \code{lambda.interp} conducts the linear interpolation of the lambdas values to obtain the coefficients at new lambda values. Note the obtained coefficients are not the exact values. \code{lamfix} fixes the largest lambda value in the lambda sequence. \code{nonzero} and \code{zeromat} organize the nonzero coefficients. Most of the aforementioned functions are modified or directly copied from the \code{gcdnet} and the \code{glmnet} packages.
}
\keyword{internal}
