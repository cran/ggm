\name{edges}
\alias{edges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Edges of a graph}
\description{
  Finds the set of edges of a graph. That is
  the set of undirected edges if the graph is undirected
  and the set of arrows if the graph is directed.
}
\usage{
edges(gmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gmat}{a square Boolean matrix, with dimnames, and edge matrix of
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
edges(UG(~ y*v*k +v*k*d+y*d))

## A DAG
edges(DAG(u~h+o+p, h~o, o~p))
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}

