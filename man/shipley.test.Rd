\name{shipley.test}
\alias{shipley.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test of all independencies implied by a DAG}
\description{
  Computes a simultaneous test of all independence relationships
  implied by a given Gaussian model  defined according to
  a directed acyclic graph, based on the sample covariance matrix.
}
\usage{
shipley.test(S, n, A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{a symmetric positive definite matrix, the sample covariance matrix.}
  \item{n}{a positive integer, the sample size.}
  \item{A}{a square Boolean matrix, of the same dimension as \code{S},
    representing the edge matrix of a DAG.}
}
}
\details{
  The test statistic is \eqn{C = -2 \sum \ln p_j} where \eqn{p_j} are the
  p-values of tests of conditional independence in the basis set
  computed by \code{basiSet(A)}. The p-values are independent
  uniform variables on \eqn{(0,1)} and the statistic has exactly a
  chi square distribution on \eqn{2k} degrees of freedom  where
  \eqn{k} is the number of elements of the basis set.
  Shipley (2002) calls this test Fisher's C test.
}
\value{
  \item{ctest}{Test statistic \eqn{C}.}
  \item{df}{Degrees of freedom.}
  \item{pvalue}{The P-value of the test, assuming a two-sided alternative.}
}
\references{Shipley, B. (2000). A new inferential test
for path models based on directed acyclic graphs. \emph{Structural
  Equation Modeling}, 7(2), 206--218.}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{basiSet}}, \code{\link{pcor.test}}}
\examples{
## A decomposable model for the mathematics marks data
data(marks)
dag <- DAG(mec ~ vec+alg, vec ~ alg, sta ~ alg+ana, ana ~ alg)
shipley.test(cov(marks), n=88, dag)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
