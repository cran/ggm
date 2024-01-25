\name{msep}
\alias{msep}
\title{The m-separation criterion
}
\description{
\code{msep} determines whether two set of nodes are m-separated
by a third set of nodes.
}
\usage{
msep(a, alpha, beta, C = c())
}
\arguments{
  \item{a}{
An adjacency matrix, or a graph that can be a \code{graphNEL} or an \code{\link{igraph}} object
or a vector of length \eqn{3e}, where \eqn{e} is the number of edges of the graph,
that is a sequence of triples (type, node1label, node2label). The type
of edge can be \code{"a"} (arrows from node1 to node2), \code{"b"} (arcs), and
\code{"l"} (lines).
}
  \item{alpha}{
A subset of the node set of \code{a}}
  \item{beta}{
Another disjoint subset of the node set of \code{a}}
  \item{C}{
A third disjoint subset of the node set of \code{a}}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A logical value. \code{TRUE}  if \code{alpha} and \code{beta} are m-separated given \code{C}.
\code{FALSE} otherwise.
}
\references{
Richardson, T.S. and Spirtes, P. (2002) Ancestral graph Markov models. \emph{Annals
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
\code{\link{dSep}}, \code{\link{MarkEqMag}}
}
\examples{
H <-matrix(c(0,0,0,0,
	         1,0,0,1,
	         0,1,0,0,
	         0,0,0,0),4,4)
msep(H,1,4, 2)
msep(H,1,4, c())
}
\keyword{graphs}
\keyword{d-separation}
\keyword{m-separation}
\keyword{mixed graph}% __ONLY ONE__ keyword per line
