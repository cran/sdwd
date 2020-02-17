\name{coef.cv.sdwd}
\alias{coef.cv.sdwd}
\title{compute coefficients from a "cv.sdwd" object}
\description{Computes coefficients at chosen values of \code{lambda} from the \code{\link{cv.sdwd}} object. 
}
\usage{
\method{coef}{cv.sdwd}(object, s=c("lambda.1se", "lambda.min"),...)
}
\arguments{
  \item{object}{A fitted \code{\link{cv.sdwd}} object, obtained by conducting the cross-validation to the sparse DWD model.}
  \item{s}{Value(s) of the L1 tuning parameter \code{lambda} for computing coefficients. Default value is \code{"lambda.1se"}, which represents the largest \code{lambda} value achieving the cross-validation error within one standard error of the minimum. An alternative value is \code{"lambda.min"}, which is the \code{lambda} incurring the least cross-validation error. \code{s} can also be numeric, being taken as the value(s) to be used.}
  \item{\dots}{Other arguments that can be passed to \code{\link{sdwd}}. } 
}

\details{This function computes the coefficients at the values of \code{lambda} suggested by the cross-validation. This function is modified based on the \code{coef.cv} function from the \code{glmnet} and the \code{gcdnet} packages.}

\value{The returned object depends on the choice of \code{s} and the \code{\dots} argument passed on to the \code{\link{sdwd}} method.}
\author{Boxiang Wang and Hui Zou\cr
Maintainer: Boxiang Wang  \email{boxiang-wang@uiowa.edu}}
\references{Wang, B. and Zou, H. (2016)
  ``Sparse Distance Weighted Discrimination", \emph{Journal of Computational and Graphical Statistics}, \bold{25}(3), 826--838.\cr
  \url{https://www.tandfonline.com/doi/full/10.1080/10618600.2015.1049700}\cr
  \cr
Yang, Y. and Zou, H. (2013)
``An Efficient Algorithm for Computing the HHSVM and Its Generalizations", 
  \emph{Journal of Computational and Graphical Statistics}, \bold{22}(2), 396--415.\cr
  \url{https://www.tandfonline.com/doi/full/10.1080/10618600.2012.680324}\cr
 \cr 
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software}, \bold{33}(1), 1--22.\cr
\url{http://www.jstatsoft.org/v33/i01/paper}
\cr}

\seealso{\code{\link{cv.sdwd}} and \code{\link{predict.cv.sdwd}} methods.}
\examples{
data(colon)
colon$x = colon$x[ , 1:100] # this example only uses the first 100 columns 
set.seed(1)
cv = cv.sdwd(colon$x, colon$y, lambda2=1, nfolds=5)
c1 = coef(cv, s="lambda.1se")
}
