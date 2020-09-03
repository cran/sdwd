\name{predict.cv.sdwd}
\alias{predict.cv.sdwd}
\title{make predictions from a "cv.sdwd" object}
\description{
This function predicts the class labels of new observations by the sparse DWD at the \code{lambda} values suggested by \code{\link{cv.sdwd}}.}
\usage{
\method{predict}{cv.sdwd}(object, newx, s=c("lambda.1se","lambda.min"),...)
}
\arguments{
    \item{object}{A fitted \code{\link{cv.sdwd}} object.}
    \item{newx}{A matrix of new values for \code{x} at which predictions are
    to be made. Must be a matrix. See documentation for \code{predict.sdwd}.}
    \item{s}{Value(s) of the L1 tuning parameter \code{lambda} for making predictions.  Default is the \code{s="lambda.1se"} saved on the \code{\link{cv.sdwd}} object. An alternative choice is \code{s="lambda.min"}. \code{s} can also be numeric, being taken as the value(s) to be used.}
    \item{\dots}{Not used. Other arguments to \code{predict}. } }
    
\details{This function uses the cross-validation results to making predictions. This function is modified based on the \code{predict.cv} function from the \code{glmnet} and the \code{gcdnet} packages.}
\value{Predicted class labels or fitted values, depending on the choice of \code{s} and the \dots argument passed on to the \code{\link{sdwd}} method.}

\author{Boxiang Wang and Hui Zou\cr
Maintainer: Boxiang Wang  \email{boxiang-wang@uiowa.edu}}
\references{
Wang, B. and Zou, H. (2016)
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
\url{https://www.jstatsoft.org/v33/i01/paper}
\cr}
 

\seealso{\code{\link{cv.sdwd}}, and \code{\link{coef.cv.sdwd}} methods.}
\examples{
data(colon)
colon$x = colon$x[ , 1:100] # this example only uses the first 100 columns 
set.seed(1)
cv = cv.sdwd(colon$x, colon$y, lambda2=1, nfolds=5)
predict(cv$sdwd.fit, newx=colon$x[2:5, ], 
  s=cv$lambda.1se, type="class")
}
