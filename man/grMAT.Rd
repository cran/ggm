\name{grMAT}
\alias{grMAT}

\title{Graph to adjacency matrix
}
\description{
\code{grMAT} generates the associated adjacency matrix to a given graph.
}
\usage{
grMAT(agr)
}
\arguments{
  \item{agr}{
A graph that can be a \code{graphNEL} or an \code{\link{igraph}} object
or a vector of length \eqn{3e}, where \eqn{e} is the number of edges of the graph,
that is a sequence of triples (type, node1label, node2label). The type
of edge can be \code{"a"} (arrows from node1 to node2), \code{"b"} (arcs), and
\code{"l"} (lines).
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A matrix that consists 4 different integers as an \eqn{ij}-element: 0 for a missing
edge between \eqn{i} and \eqn{j}, 1 for an arrow from \eqn{i} to \eqn{j}, 10 for a full line between
\eqn{i} and \eqn{j}, and 100 for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers are
added to be associated with multiple edges of different types. The matrix is
symmetric w.r.t full lines and bi-directed arrows.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Kayvan Sadeghi
}
%\note{
%%  ~~further notes~~
%}


%\seealso{
%\code{\link{adjGR}}
%}
\examples{
## Generating the adjacency matrix from a vector
exvec <-c ('b',1,2,'b',1,14,'a',9,8,'l',9,11,'a',10,8,
           'a',11,2,'a',11,10,'a',12,1,'b',12,14,'a',13,10,'a',13,12)
grMAT(exvec)
}
\keyword{graphs}
\keyword{adjacency matrix}
\keyword{mixed graph}
\keyword{vector}

