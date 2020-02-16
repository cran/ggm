\name{parcor}
\alias{parcor}
\title{Partial correlations}
\description{
  Finds the matrix of the partial correlations between pairs
  of variables given the rest.
}
\usage{
parcor(S)
}
\arguments{
  \item{S}{a symmetric positive definite matrix, representing a
    covariance matrix.}
}
\details{
  The algorithm computes \eqn{- \sigma^{rs}/(\sigma^{rr}
    \sigma^{ss})^{1/2}} where the \eqn{\sigma^{rs}} are concentrations,
    i.e. elements of the inverse covariance matrix.  
}
\value{
  A symmetric matrix with ones along the diagonal and in position
  \eqn{(r,s)} the partial correlation between variables \eqn{r}
  and \eqn{s} given all the remaining variables.}
\references{Cox, D. R. \& Wermuth, N. (1996). \emph{Multivariate
    dependencies}. London: Chapman \& Hall.}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{var}}, \code{\link{cor}}, \code{\link{correlations}}}
\examples{
### Partial correlations for the mathematics marks data
data(marks)
S <- var(marks)
parcor(S)
}
\keyword{array}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}

