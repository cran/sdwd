\name{sdwd}
\alias{sdwd}
\title{fit the sparse DWD}
\description{Fits the sparse distance weighted discrimination (DWD) model with imposing L1, elastic-net, or adaptive elastic-net penalties. The solution path is computed at a grid of values of tuning parameter \code{lambda}. This function is modified based on the \code{glmnet} and the \code{gcdnet} packages.}
\usage{
sdwd(x, y, nlambda=100, 
     lambda.factor=ifelse(nobs < nvars, 0.01, 1e-04), 
     lambda=NULL, lambda2=0, pf=rep(1, nvars), 
     pf2=rep(1, nvars), exclude, dfmax=nvars + 1, 
     pmax=min(dfmax * 1.2, nvars), standardize=TRUE, 
     eps=1e-8, maxit=1e6, strong=T)
}
\arguments{
		\item{x}{A matrix with \eqn{N} rows and \eqn{p} columns for predictors.}
		\item{y}{A vector of length \eqn{p} for binary responses. The element of \code{y} is either -1 or 1.}
		\item{nlambda}{The number of \code{lambda} values, i.e., length of the \code{lambda} sequence. Default is 100.}
		\item{lambda.factor}{The ratio of the smallest to the largest \code{lambda} in the sequence: \code{lambda.factor} = \code{min(lambda)} / \code{max(lambda)}.  \code{max(lambda)} is the least \code{lambda} to make all coefficients to be zero. The default value of \code{lambda.factor} is 0.0001 if \eqn{N >= p} or 0.01 if \eqn{N < p}. Takes no effect when user specifies a \code{lambda} sequence.} 
		\item{lambda}{An optional user supplied \code{lambda} sequence. If \code{lambda = NULL} (default), the program computes its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.factor}; otherwise, the program uses the user-specified one. Since the program will automatically sort user-defined \code{lambda} sequence in decreasing order, it is better to supply a decreasing sequence.}
		\item{lambda2}{The L2 tuning parameter \eqn{\lambda_2}{lambda2}.}
		\item{pf}{A vector of length \eqn{p}{p} representing the L1 penalty weights to each coefficient of \eqn{\beta}{beta} for adaptive L1 or adaptive elastic net. \code{pf} can be 0 for some predictor(s), leading to including the predictor(s) all the time. One suggested choice of \code{pf} is \eqn{{\beta + 1/n}^{-1}}{1/(beta+1/n)}, where \eqn{n} is the sample size and \eqn{\beta}{beta} is the coefficents obtained by L1 DWD or enet DWD. Default is 1 for all predictors (and infinity if some predictors are listed in \code{exclude}).}
		\item{pf2}{A vector of length \eqn{p}{p} for L2 penalty factor for adaptive L1 or adaptive elastic net. To allow different L2 shrinkage, user can set \code{pf2} to be different L2 penalty weights for each coefficient of \eqn{\beta}{beta}. \code{pf2} can be 0 for some variables, indicating no L2 shrinkage. Default is 1 for all predictors.}
		\item{exclude}{Whether to exclude some predictors from the model. This is equivalent to adopting an infinite penalty factor when excluding some predictor. Default is none. }
		\item{dfmax}{Restricts at most how many predictors can be incorporated in the model. Default is \eqn{p+1}. This restriction is helpful when \eqn{p} is large, provided that a partial path is acceptable. }
		\item{pmax}{Restricts the maximum number of variables ever to be nonzero; e.g, once some \eqn{\beta} enters the model, it counts once. The count will not change when the \eqn{\beta} exits or re-enters the model. Default is \code{min(dfmax*1.2,p)}.}
		\item{standardize}{Whether to standardize the data. If \code{TRUE}, \code{\link{sdwd}} normalizes the predictors such that each column has sum squares\eqn{\sum^N_{i=1}x_{ij}^2/N=1}{<Xj,Xj>/N} of one. Note that x is always centered (i.e. \eqn{\sum^N_{i=1}x_{ij}=0}{sum(Xj)=0}) no matter \code{standardize} is \code{TRUE} or \code{FALSE}. \code{\link{sdwd}} always returns coefficient \code{beta} on the original scale.  Default value is \code{TRUE}.}
		\item{eps}{The algorithm stops when (i.e. \eqn{4\max_j(\beta_j^{new}-\beta_j^{old})^2}{4*max(j)(beta_new[j]-beta_old[j])^2} is less than \code{eps}, where \eqn{j=0,\ldots, p}. Defaults value is \code{1e-8}.}
		\item{maxit}{Restricts how many outer-loop iterations are allowed. Default is 1e6. Consider increasing \code{maxit} when the algorithm does not converge.}
		\item{strong}{If \code{TRUE}, adopts the strong rule to accelerate the algorithm.}
}

\details{
The \code{\link{sdwd}} minimizes the sparse penalized DWD loss function, 
\deqn{L(y, X, \beta)/N + \lambda_1||\beta||_1 + 0.5\lambda_2||\beta||_2^2,}{L(y, X, beta))/N + lambda1 * ||beta||_1^1 + 0.5 * lambda2 * ||beta||_2^2,}
where \eqn{L(u)=1-u}{L(u) = 1 - u} if \eqn{u \le 1/2}{u <= 0.5}, \eqn{1/(4u)}{1 / (4*u)} if \eqn{u > 1/2}{u > 0.5} is the DWD loss. The value of \code{lambda2} is user-specified.

To use the L1 penalty (lasso), set \code{lambda2=0}. To use the elastic net, set \code{lambda2} as nonzero. To use the adaptive L1, set \code{lambda2=0} and specify \code{pf} and \code{pf2}. To use the adaptive elastic net, set \code{lambda2} as nonzero and specify \code{pf} and \code{pf2} as well.  

When the algorithm do not converge or run slow, consider increasing \code{eps}, decreasing
\code{nlambda}, or increasing \code{lambda.factor} before increasing
\code{maxit}.
}

\value{
An object with S3 class \code{\link{sdwd}}.
		\item{b0}{A vector of length \code{length(lambda)} representing the intercept at each \code{lambda} value.}
		\item{beta}{A matrix of dimension \code{p*length(lambda)} representing the coefficients at each \code{lambda} value. The matrix is stored as a sparse matrix  (\code{Matrix} package). To convert it into normal type matrix use \code{as.matrix()}.}
		\item{df}{The number of nonzero coefficients at each \code{lambda}.}
		\item{dim}{The dimension of coefficient matrix, i.e., \code{p*length(lambda)}.}
		\item{lambda}{The \code{lambda} sequence that was actually used.}
        \item{npasses}{Total number of iterations for all lambda values. }
		\item{jerr}{Warnings and errors; 0 if no error.}
    \item{call}{The call which produced this object.}
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
\url{http://www.jstatsoft.org/v33/i01/paper}\cr
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

\seealso{\code{print.sdwd}, \code{predict.sdwd}, \code{coef.sdwd}, \code{plot.sdwd}, and \code{cv.sdwd}.}
\examples{
# load the data
data(colon)
# fit the elastic-net penalized DWD with lambda2=1
fit = sdwd(colon$x, colon$y, lambda2=1)
print(fit)
# coefficients at some lambda value
c1 = coef(fit, s=0.005)
# make predictions
predict(fit, newx=colon$x[1:10, ], s=c(0.01, 0.005))

}
\keyword{models}
\keyword{regression}
