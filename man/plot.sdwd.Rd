\name{plot.sdwd}
\alias{plot.sdwd}
\title{plot coefficients for the sparse DWD}
\description{
Plots the solution paths for a fitted \code{\link{sdwd}} object.}
\usage{
\method{plot}{sdwd}(x, xvar=c("norm", "lambda"), color=FALSE, label=FALSE, ...)
}
\arguments{
    \item{x}{A fitted \code{\link{sdwd}} model.}
    \item{xvar}{Specifies the X-axis. If \code{xvar == "norm"}, plots against the  L1-norm of the coefficients; if \code{xvar == "lambda"} against the log-lambda  sequence.}
    \item{color}{If \code{TRUE}, plots the curves with rainbow colors; otherwise, with gray colors (default).}
    \item{label}{If \code{TRUE}, labels the curves with variable sequence numbers. Default is \code{FALSE}. }
    \item{\dots}{Other graphical parameters to \code{plot}.}
}
\details{
Plots the solution paths as a coefficient profile plot. This function is modified based on the \code{plot} function from the \code{gcdnet} and the \code{glmnet} packages.
}
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
}
 
 \seealso{\code{print.sdwd}, \code{predict.sdwd}, \code{coef.sdwd}, \code{plot.sdwd}, and \code{cv.sdwd}.}
 
\examples{
data(colon)
fit = sdwd(colon$x, colon$y)
par(mfrow=c(1,3))
# plots against the L1-norm of the coefficients
plot(fit) 
# plots against the log-lambda sequence
plot(fit, xvar="lambda", label=TRUE)
# plots with colors
plot(fit, color=TRUE)
}
