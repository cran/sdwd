\name{coef.sdwd}
\alias{coef.sdwd}
\title{compute coefficients for the sparse DWD}
\description{
Computes the coefficients or returns the indices of nonzero coefficients at chosen values of \code{lambda} from a fitted \code{\link{sdwd}} object.
}
\usage{
\method{coef}{sdwd}(object, s=NULL, type=c("coefficients","nonzero"), ...)
}
\arguments{
		\item{object}{A fitted \code{\link{sdwd}} object.}
		\item{s}{Value(s) of the L1 tuning parameter \code{lambda} for computing coefficients. Default is the entire \code{lambda} sequence obtained by \code{\link{sdwd}}.}
		\item{type}{\code{"coefficients"} or \code{"nonzero"}? \code{"coefficients"} computes the coefficients at given values for \code{s}; \code{"nonzero"} returns a list of the indices of the nonzero coefficients for each value of \code{s}. Default is \code{"coefficients"}.}
		\item{\dots}{Not used. Other arguments to predict. } 
}

\details{
\code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda sequence used for fitting the model, the \code{coef} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of coefficients from both left and right \code{lambda} indices. This function is modified based on the \code{coef} function from the \code{gcdnet} and the \code{glmnet} packages.
}
\value{Either the coefficients at the requested values of \code{lambda}, or a list of the indices of the nonzero coefficients for each \code{lambda}.}


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
\cr

 }
 

\seealso{\code{\link{predict.sdwd}}}
\examples{
data(colon)
fit = sdwd(colon$x, colon$y, lambda2=1)
c1 = coef(fit, type="coef",s=c(0.1, 0.005))
c2 = coef(fit, type="nonzero")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
