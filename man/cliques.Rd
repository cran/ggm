\name{cliques}
\alias{cliques}
\title{Cliques of an undirected graph}
\description{
Finds the cliques of an undirected graph.
}
\usage{
cliques(gmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gmat}{a square Boolean matrix with dimnames, representing the
  edge matrix of an undirected graph.}
}
\details{
  The cliques of a graph are the subsets of nodes which induce a
  maximally complete subgraph. The algorithm is NP-complete in general.
}
\value{
  a list of vectors of integers representing  subsets of the set of nodes. 
}
\references{
  Lauritzen, S. (1996). \emph{Graphical models}. Oxford:
  Clarendon Press.}
\author{The code was written by somebody at  Lancaster University.}
\seealso{\code{\link{UG}},\code{\link{bd}}}
\examples{
## Two cliques
u <- UG(~ a*b*c + c*d*e*g)
u
cliques(u)

## the same (note that the formula is not of a graphical model)
u <- UG(~ a*b+b*c+c*a + c*d*e*g)
u
cliques(u)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}

