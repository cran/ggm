\name{triDec}
\alias{triDec}
\title{Triangular decomposition of a covariance matrix}
\description{
  Decomposes a symmetric positive definite matrix with a variant of the
  Cholesky decomposition. 
}
\usage{
triDec(Sigma)
}
\arguments{
  \item{Sigma}{a symmetric positive definite matrix.}
}
\details{
  Any symmetric positive definite matrix \eqn{\Sigma}{Sigma}
  can be decomposed as
  \eqn{\Sigma = B \Delta B^T}{Sigma = B \%*\% Delta \%*\% t(B)}
  where \eqn{B} is upper triangular with ones
  along the main diagonal and \eqn{\Delta}{Delta} is diagonal. If
  \eqn{\Sigma}{Sigma} is a covariance
  matrix, the concentration matrix is \eqn{\Sigma^{-1} = A^T \Delta^{-1}
    A} where  \eqn{A = B^{-1}} is the matrix of the regression coefficients
  (with the sign changed) of a system of linear recursive regression
  equations  with independent residuals. In the equations each variable 
  \eqn{i} is regressed on the variables \eqn{i+1, \dots, d}.
  The elements on the diagonal of \eqn{\Delta} are the partial variances.
}
\value{
  \item{A}{a square upper triangular matrix of the same order as
    \code{Sigma} with ones on the diagonal.}
  \item{B}{the inverse of \code{A}, another triangular matrix with unit diagonal.}
  \item{Delta}{a vector containing the diagonal values of \eqn{\Delta}.}
}
\references{Cox, D. R. \& Wermuth, N. (1996). \emph{Multivariate
  dependencies}. London: Chapman \& Hall.}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{chol}}}
\examples{
## Triangular decomposition of a covariance matrix
B <- matrix(c(1,  -2, 0, 1,
              0,   1, 0, 1,
              0,   0, 1, 0,
              0,   0, 0, 1), 4, 4, byrow=TRUE)
B
D <- diag(c(3, 1, 2, 1))
S <- B \%*\% D \%*\% t(B)
triDec(S)
solve(B)
}
\keyword{array}
\keyword{algebra}
\keyword{models}
\keyword{multivariate}
