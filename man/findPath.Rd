\name{findPath}
\alias{findPath}
\title{Finding paths}
\description{
Finds one path between two nodes of a graph.
}
\usage{
findPath(amat, st, en, path = c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{a square Boolean matrix with dimnames, the adjacency matrix of a graph.}
  \item{st}{an integer, the starting node.}
  \item{en}{an integer, the ending node.}
  \item{path}{a vector of integers, used in recursive calls. At the
    beginning is \code{NULL}. It should not be modified by the user.}
}
\value{
  a vector of integers, the sequence of nodes of a path,
  starting from \code{st} to \code{en}. In some graphs
  (spanning trees) there is only one path between two nodes.
}
\note{
  This function is not intended to be directly called by the user.
}
\references{Python Softftware Foundation (2003). Python Patterns ---   Implementing Graphs. \url{https://www.python.org/doc/essays/graphs/}.}
\author{Giovanni M. Marchetti, translating the original \pkg{Python}
  code (see references).}
\seealso{\code{\link{fundCycles}}}
\examples{
## A (single) path on a spanning tree
findPath(bfsearch(UG(~ a*b*c + b*d + d*e+ e*c))$tree, st=1, en=5)
}
\keyword{graphs}
