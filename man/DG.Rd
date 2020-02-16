\name{DG}
\alias{DG}
\title{Directed graphs}
\description{
  Defines the adjacency of a directed graph.
  }
\usage{
DG(\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{a sequence of model formulae}
}
\details{
  The directed graph is defined by a sequence of models formulae.
  For each formula the response defines a node of the graph and
  its parents. The graph contains no loops.  
}
\value{the adjacency matrix of the directed graph, i.e., 
  a square Boolean matrix of order equal to the number of nodes of the
  graph and a one in position \eqn{(i,j)} if there is an arrow from
  \eqn{i} to \eqn{j} and zero otherwise.  The dimnames of the adjacency
  matrix are the labels for the nodes of the graph. 
}
\references{
  Lauritzen, S. (1996). \emph{Graphical models}. Oxford:
  Clarendon Press.
}
\author{G. M. Marchetti}
\seealso{\code{\link{DAG}}, \code{\link{UG}}}
\examples{
## A DAG
DG(y ~ x, x ~ z, z ~ u)

## A cyclic directed graph
DG(y ~ x, x ~ z, z ~ y)

## A graph with two arrows between two nodes
DG(y ~ x, x ~ y)

## There can be isolated nodes
DG(y ~ x, x ~ x)
}
\keyword{graphs}   
\keyword{directed graph}
\keyword{models}
\keyword{multivariate}
