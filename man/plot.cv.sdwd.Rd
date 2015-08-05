\name{plot.cv.sdwd}
\alias{plot.cv.sdwd}
\title{plot the cross-validation curve of the sparse DWD}
\description{
Plots the cross-validation curve against a function of \code{lambda} values. The function also provides the upper and lower standard deviation curves.}
\usage{
\method{plot}{cv.sdwd}(x, sign.lambda, ...)
}
\arguments{
	\item{x}{A fitted \code{\link{cv.sdwd}} object.}
	\item{sign.lambda}{Whether to plot against \code{log(lambda)} (default) or its negative if \code{sign.lambda=-1}.}
	\item{\dots}{Other graphical parameters to \code{plot}.}
}
\details{This function depicts the cross-validation curves. This function is modified based on the \code{plot.cv} function from the \code{glmnet} and the \code{gcdnet} packages.}

\author{Boxiang Wang and Hui Zou\cr
Maintainer: Boxiang Wang  \email{boxiang@umn.edu}}
\references{Wang, B. and Zou, H. (2015)
  ``Sparse Distance Weighted Discrimination", \emph{Journal of Computational and Graphical Statistics}, forthcoming. \cr
  \url{http://arxiv.org/abs/1501.06066}\cr
  \cr
Yang, Y. and Zou, H. (2013)
``An Efficient Algorithm for Computing the HHSVM and Its Generalizations", 
  \emph{Journal of Computational and Graphical Statistics}, \bold{22}(2), 396--415\cr
  \url{http://users.stat.umn.edu/~yiyang/resources/papers/JCGS_gcdnet.pdf}\cr
 \cr 
  Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software}, \bold{33}(1), 1--22\cr
\url{http://www.jstatsoft.org/v33/i01/paper}
\cr}
 
\seealso{\code{\link{cv.sdwd}}.}
\examples{
data(colon)
colon$x = colon$x[ ,1:100] # this example only uses the first 100 columns 
set.seed(1)
cv = cv.sdwd(colon$x, colon$y, lambda2=1, nfolds=5)
plot(cv)
}
\keyword{models}
\keyword{regression}
