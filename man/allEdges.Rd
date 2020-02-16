\name{allEdges}
\alias{allEdges}
\title{All edges of a graph}
\description{
  Finds the set of edges of a graph. That is
  the set of undirected edges if the graph is undirected
  and the set of arrows if the graph is directed.
}
\usage{
allEdges(amat)
}
\arguments{
  \item{amat}{a square Boolean matrix, with dimnames, the adjacency matrix of
    a graph.}
}
\value{
  a matrix with two columns. Each row of the matrix is a pair of
  indices indicating an edge of the graph. If the graph is undirected,
  then only one of the pairs \eqn{(i,j), (j,i)} is reported.
}
\author{Giovanni M. Marchetti} 
\seealso{\code{\link{cycleMatrix}}}
\examples{
## A UG graph
allEdges(UG(~ y*v*k +v*k*d+y*d))

## A DAG
allEdges(DAG(u~h+o+p, h~o, o~p))
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}

