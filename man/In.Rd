\name{In}
\alias{In}
\title{Indicator matrix}
\description{
  Finds the indicator matrix of the zeros of a matrix.
}
\usage{
In(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{a matrix.}
}
\details{
The indicator matrix is a matrix of zeros and ones which has a zero
element iff the corresponding element of \code{A} is (exactly) zero.
}
\value{a matrix of the same dimensions as \code{A}.
}
\references{
  Wermuth, N. & Cox, D.R. (2004). Joint response graphs
  and separation induced by triangular systems. \emph{J.R. Statist. Soc. B},
  66, Part 3, 687-717.
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{DAG}}, \code{\link{inducedCovGraph}}, \code{\link{inducedConGraph}}}
\examples{
## A simple way to find the overall induced concentration graph
## The DAG on p. 198 of Cox & Wermuth (1996)
amat <- DAG(y1 ~ y2 + y3, y3 ~ y5, y4 ~ y5)
A <- edgematrix(amat)
In(crossprod(A))
}
\keyword{array}
\keyword{algebra}
\keyword{graphs}
\keyword{multivariate}
