\name{is.acyclic}
\alias{is.acyclic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graph queries}
\description{
  Checks if a given graph is acyclic.
}
\usage{
is.acyclic(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{a square Boolean matrix with dimnames, the edge matrix of a graph.}
}
\value{
  a logical value, \code{TRUE} if the graph is acyclic and \code{FALSE}
  otherwise.
}
\references{
  Aho, A.V., Hopcrtoft, J.E. \& Ullman, J.D. (1983). \emph{Data structures
    and algorithms.} Reading: Addison-Wesley.
}
\author{Giovanni M. Marchetti} 
\seealso{\code{\link{ancGraph}}}
\examples{
## Force a cyclic graph
d <- DAG(y ~ x, x ~ z, z~ y, test=FALSE)
## Test acyclicity
is.acyclic(d)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
