\name{Max}
\alias{Max}
\title{
Maximisation for graphs
}
\description{
\code{Max} generates a maximal graph that induces the same
independence model from a non-maximal graph.
}
\usage{
Max(amat)
}
\arguments{
  \item{amat}{
An adjacency matrix, or a graph that can be a \code{graphNEL} or an \code{\link{igraph}} object
or a vector of length \eqn{3e}, where \eqn{e} is the number of edges of the graph,
that is a sequence of triples (type, node1label, node2label). The type
of edge can be \code{"a"} (arrows from node1 to node2), \code{"b"} (arcs), and
\code{"l"} (lines).
}
}
\details{
\code{Max} looks for non-adjacent pais of nodes that are connected by
primitive inducing paths, and connect such pairs by an appropriate edge.
}
\value{
A matrix that consists 4 different integers as an \eqn{ij}-element: 0 for a missing
edge between \eqn{i} and \eqn{j}, 1 for an arrow from \eqn{i} to \eqn{j}, 10 for a full line between
\eqn{i} and \eqn{j}, and 100 for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers are
added to be associated with multiple edges of different types. The matrix is
symmetric w.r.t full lines and bi-directed arrows.
}
\references{
Richardson, T.S. and Spirtes, P. (2002). Ancestral graph Markov models. \emph{Annals
of Statistics}, 30(4), 962-1030.

Sadeghi, K. and Lauritzen, S.L. (2014). Markov properties for loopless mixed graphs. \emph{Bernoulli} 20(2), 676-696.
}
\author{
Kayvan Sadeghi
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{MAG}}, \code{\link{MRG}}, \code{\link{msep}}, \code{\link{MSG}}
}
\examples{
H <- matrix(c(  0,100,  1,  0,
	          100,  0,100,  0,
	            0,100,  0,100,
	            0,  1,100,  0), 4, 4)
Max(H)
}
\keyword{graphs}
\keyword{loopless mixed graph}
\keyword{m-separation}
\keyword{maximality}

