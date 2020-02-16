\name{RepMarUG}
\alias{RepMarUG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Representational Markov equivalence to undirected graphs.
}
\description{
\code{RepMarUG} determines whether a given maximal ancestral graph can 
be Markov equivalent to an undirected graph, and if that is the case, it finds
an undirected graph that is Markov equivalent to the given graph.}
\usage{
RepMarUG(amat)
}
%- maybe also 'usage' for other objects documented here.
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
\code{RepMarBG} looks for presence of an unshielded collider 
V-configuration in graph.}
\value{
  A list with two components: \code{verify} and
  \code{amat}. \code{verify} is a logical value, \code{TRUE} if there is
  a representational Markov equivalence and \code{FALSE} otherwise.
  \code{amat} is either \code{NA} if \code{verify == FALSE} or
  the adjacency matrix of the generated graph, if
  \code{verify == TRUE}. In this case it consists of 4 different
  integers as an \eqn{ij}-element: 0 for a missing 
edge between \eqn{i} and \eqn{j}, 1 for an arrow from \eqn{i} to \eqn{j}, 10 for a full line between 
\eqn{i} and \eqn{j}, and 100 for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers are 
added to be associated with multiple edges of different types. The matrix is 
symmetric w.r.t full lines and bi-directed arrows.
}
\references{
Sadeghi, K. (2011). Markov equivalences for subclasses of loopless mixed graphs.
\emph{Submitted}, 2011. }
\author{
Kayvan Sadeghi
}
%\note{
%%  ~~further notes~~
%}


\seealso{
\code{\link{MarkEqMag}}, \code{\link{MarkEqRcg}}, \code{\link{RepMarBG}},
\code{\link{RepMarDAG}}
}
\examples{
H<-matrix(c(0,10,0,0,10,0,0,0,0,1,0,100,0,0,100,0),4,4)
RepMarUG(H)
}
\keyword{graphs}
\keyword{bidirected graph}
\keyword{Markov equivalence}
\keyword{maximal ancestral graph}
\keyword{representational Markov equivalence}
