\name{topSort}
\alias{topSort}
\alias{topOrder}
\title{Topological sort}
\description{
  \code{topOrder} returns the topological order of a directed
  acyclic graph (parents, before children). \code{topSort} permutates
  the adjacency matrix according to the topological order.
}
\usage{
topSort(amat)
topOrder(amat)
}
\arguments{
  \item{amat}{a square Boolean matrix with dimnames, representing the
    adjacency matrix of a directed acyclic graph.}
}
\details{
  The topological order needs not to be unique.
  After the permutation the adjacency matrix of the graph is
  upper triangular. The function is a translation of
  the Matlab function \code{topological_sort} in Toolbox
  \pkg{BNT} written by Kevin P. Murphy.
}
\value{
  \code{topOrder(amat)} returns a vector of integers representing
  the permutation of the nodes. \code{topSort(amat)} returns
  the adjacency matrix with rows and columns permutated.
}
\note{
  The order of the nodes defined by \code{DAG} is that of their first
  appearance in the model formulae (from left to right).
}
\references{
  Aho, A.V., Hopcrtoft, J.E. & Ullman, J.D. (1983). \emph{Data structures
    and algorithms.} Reading: Addison-Wesley.

  Lauritzen, S. (1996). \emph{Graphical models}. Oxford:
  Clarendon Press.
}
\author{Kevin P. Murphy, Giovanni M. Marchetti}
\seealso{\code{\link{DAG}}, \code{\link{isAcyclic}}}
\examples{
## A simple example
dag <- DAG(a ~ b, c ~ a + b, d ~ c + b)
dag
topOrder(dag)
topSort(dag)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
