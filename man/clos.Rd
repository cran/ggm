\name{clos}
\alias{clos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graph operations}
\description{
Transforms an edge matrix of a graph.
}
\usage{
clos(M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{a square matrix with dimnames representing the edge matrix of
  an undirected graph}
}
\details{
  The function \code{clos} is not intended to be used directly
  from the user. It is used by functions that compute the induced
  covariance or concentration graphs. The function \code{clos}
  is explained by Cox \& Wermuth (2003).
}
\value{
  a Boolean matrix of the same dimensions of \code{M}.}
\references{
  Wermuth, N. \& Cox, D.R. (2003). Joint response graphs
  and separation induced by triangular systems. Submitted and available at
  \url{http://psystat.sowi.uni-mainz.de}.
  }
\author{Giovanni M. Marchetti} 
\seealso{\code{\link{DAG}}, \code{\link{inducedCovGraph}}, \code{\link{inducedConGraph}}}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}

