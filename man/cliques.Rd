\name{cliques}
\alias{cliques}
\title{Cliques of an undirected graph}
\description{
Finds the cliques of an undirected graph.
}
\usage{
cliques(amat)
}
\arguments{
  \item{amat}{a square Boolean matrix with dimnames, representing the
  adjacency matrix of an undirected graph.}
}
\details{
  The cliques of a graph are the subsets of nodes which induce a
  maximally complete subgraph. Determining the cliques of a graph is
  NP-hard in general. 
}
\value{
  a list of vectors of nodes. 
}
\references{
  Lauritzen, S. (1996). \emph{Graphical models}. Oxford:
  Clarendon Press.}
\author{Mathias Drton}
\seealso{\code{\link{UG}}, \code{\link{bd}}}
\examples{
u <- UG(~ a*b*c + c*d*e*g)
u
cliques(u)

graph22 <- UG(~x1*x2*x3+x3*x4*x5*x6+x5*x7*x8*x9*x10+x9*x11*x12*x13*x14+
x12*x15+x15*x16*x17*x18+x17*x1*x18+x1*x5*x7*x19*x20+x6*x20*x21+x22)

cliques(graph22)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}

