\name{rnormDag}
\alias{rnormDag}
\title{Random sample from a decomposable Gaussian model}
\description{
  Generates a sample from a mean centered  multivariate normal
  distribution whose covariance matrix has a given triangular
  decomposition.
}
\usage{
rnormDag(n, A, Delta)
}
\arguments{
  \item{n}{an integer > 0, the sample size.}
  \item{A}{a square, upper triangular matrix with ones along the
    diagonal. It defines, together with \code{Delta},
    the concentration matrix (and also the covariance matrix)
    of the multivariate normal. The order of \code{A} is the
    number of components of the normal.}
  \item{Delta}{a numeric vector of length equal to the number of columns
    of \code{A}.}
}
\details{
  The value in position \eqn{(i,j)} of  \code{A} (with \eqn{i < j}) is
  a regression coefficient (with  sign changed) in the regression of
  variable \eqn{i} on variables \eqn{i+1, \dots, d}.

  The value in position \eqn{i} of \code{Delta} is the residual
  variance in the above regression.
  
}
\value{a matrix with \code{n} rows and \code{nrow(A)} columns,
  a sample from a multivariate normal distribution with mean zero
  and covariance matrix
  \code{S = solve(A) \%*\% diag(Delta) \%*\% t(solve(A))}.
}
\references{Cox, D. R. \& Wermuth, N. (1996). \emph{Multivariate
  dependencies}. London: Chapman \& Hall.}
\author{Giovanni M. Marchetti}

\seealso{\code{\link{triDec}}, \code{\link{fitDag}}}
\examples{
## Generate a sample of 100 observation from a multivariate normal
## The matrix of the path coefficients 
A <- matrix(
c(1, -2, -3,  0, 0,  0,  0,
  0,  1,  0, -4, 0,  0,  0,
  0,  0,  1,  2, 0,  0,  0,
  0,  0,  0,  1, 1, -5,  0,
  0,  0,  0,  0, 1,  0,  3,
  0,  0,  0,  0, 0,  1, -4,
  0,  0,  0,  0, 0,  0,  1), 7, 7, byrow=TRUE)
D <- rep(1, 7)
X <- rnormDag(100, A, D)

## The true covariance matrix
solve(A) \%*\% diag(D) \%*\% t(solve(A))

## Triangular decomposition of the sample covariance matrix
triDec(cov(X))$A
}
\keyword{distribution}
\keyword{multivariate}
