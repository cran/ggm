\name{MarkEqRcg}
\alias{MarkEqRcg}
\title{
Markov equivalence for regression chain graphs.
}
\description{
\code{MarkEqMag} determines whether two RCGs (or subclasses of RCGs) are 
Markov equivalent.
}
\usage{
MarkEqRcg(amat, bmat)
}
\arguments{
  \item{amat}{
An adjacency matrix of an RCG or a graph that can be a \code{graphNEL} or an \code{\link{igraph}} object 
or a vector of length \eqn{3e}, where \eqn{e} is the number of edges of the graph, 
that is a sequence of triples (type, node1label, node2label). The type 
of edge can be \code{"a"} (arrows from node1 to node2), \code{"b"} (arcs), and 
\code{"l"} (lines).
}
  \item{bmat}{
The same as \code{amat}.
}
}
\details{
The function checks whether the two graphs have the same skeleton and 
unshielded colliders.
}
\value{
"Markov Equivalent" or "NOT Markov Equivalent". 
}
\references{
Wermuth, N. and Sadeghi, K. (2011). Sequences of regressions and their independences.
\emph{TEST}, To appear.
\url{http://arxiv.org/abs/1103.2523}.
}
\author{
Kayvan Sadeghi
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{MarkEqMag}}, \code{\link{msep}}
}
\examples{
H1<-matrix(c(0,100,0,0,0,100,0,100,0,0,0,100,0,0,0,1,0,0,0,100,0,0,1,100,0),5,5)
H2<-matrix(c(0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,100,0,0,1,100,0),5,5)
H3<-matrix(c(0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0),5,5)
#MarkEqRcg(H1,H2)
#MarkEqRcg(H1,H3)
#MarkEqRcg(H2,H3)
}
\keyword{graphs}
\keyword{bidirected graph}
\keyword{directed acyclic graph}
\keyword{Markov equivalence}
\keyword{regression chain graph}
\keyword{undirected graph}  
\keyword{multivariate}