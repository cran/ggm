\name{correlations}
\alias{correlations}
\title{Marginal and partial correlations}
\description{
  Computes a correlation matrix with ones along the
  diagonal, marginal correlations in the lower triangle and
  partial correlations given all remaining variables in the
  upper triangle.
}
\usage{
correlations(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a square symmetric matrix, a covariance matrix, or a
  data.frame for n observations and p variables.}
}
\value{
  a square correlation matrix with marginal correlations (lower
  triangle) and partial correlations (upper triangle).
}
\references{Cox, D. R. \& Wermuth, N. (1996). \emph{Multivariate
  dependencies}. London: Chapman \& Hall.}
\author{Giovanni M. Marchetti} 
\seealso{\code{\link{parcor}},  \code{\link{cor}}}
\examples{
## See Table 6.1 in Cox & Wermuth (1996)
data(glucose)
correlations(glucose)
}
\keyword{array}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}

