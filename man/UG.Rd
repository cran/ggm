\name{UG}
\alias{UG}
\title{Defining an undirected graph (UG)}
\description{
  A simple way to define an undirected graph by means of a single
  model formula.
}
\usage{
UG(f)
}
\arguments{
  \item{f}{a single model formula without response}
}
\details{
  The undirected graph \eqn{G = (V, E)} is defined by a set of nodes
  \eqn{V} and a set of pairs \eqn{E}. The set of pairs is defined by
  the set of interactions in the formula. Interactions 
  define complete subgraphs (not necessarily maximal) of the UG.
  The best way is to specify interactions that match the cliques
  of the undirected graph. This is the standard way to define
  graphical models for contingency tables. Remember that some
  hierarchical models are not graphical, but they imply the same graph.

  The function returns the edge matrix of the graph, i.e. 
  a square Boolean matrix of order equal to the number of nodes of the
  graph and a one in position \eqn{(i,j)} if there is an arrow from
  \eqn{j} to \eqn{i} and zero otherwise. By default this matrix
  has ones along the main diagonal. For UGs this matrix is symmetric.
  The dimnames of the edge matrix are the nodes of the UG.
}
\value{
  a Boolean matrix with dimnames, 
  the adjacency matrix of the undirected graph.
}
\references{Lauritzen, S. (1996). \emph{Graphical models}. Oxford:
Clarendon Press.}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{fitConGraph}}, \code{\link{fitCovGraph}}, \code{\link{DAG}}}
\examples{
## X independent of Y given Z
UG(~ X*Z + Y*Z)

# The saturated model
UG(~ X*Y*Z)

## The model without three-way interactions has the same graph
UG(~ X*Y + Y*Z + Z*X)
UG(~ (X + Y + Z)^2)

## Butterfly model defined from the cliques
UG(~ mec*vec*alg + alg*ana*sta)

## Some isolated nodes
UG(~x*y*z + a + b) 
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
