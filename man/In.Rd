\name{In}
\alias{In}
%- Also NEED an '\alias' for EACH other topic documented here.
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
  Wermuth, N. \& Cox, D.R. (2003). Joint response graphs
  and separation induced by triangular systems. Submitted and available at
  \url{http://psystat.sowi.uni-mainz.de}.
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{DAG}}, \code{\link{inducedCovGraph}}, \code{\link{inducedConGraph}}}
\examples{
## A simple way to find the overall induced coincentration graph
## The DAG on p. 198 of Cox & Wermuth (1996)
A <- DAG(y1 ~ y2 + y3, y3 ~ y5, y4 ~ y5)
In(crossprod(A))
}
\keyword{array}
\keyword{algebra}
\keyword{graphs}
\keyword{multivariate}
