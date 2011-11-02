\name{MarkEqMag}
\alias{MarkEqMag}
\title{
Markov equivalence of maximal ancestral graphs 
}
\description{
\code{MarkEqMag} determines whether two MAGs are Markov equivalent.
}
\usage{
MarkEqMag(amat, bmat)
}

\arguments{
  \item{amat}{
An adjacency matrix of a MAG, or a graph that can be a \code{graphNEL} or an \code{\link{igraph}} object 
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
 colliders with order.
}
\value{
"Markov Equivalent" or "NOT Markov Equivalent". 
}
\references{
Ali, R.A.,  Richardson, T.S. and Spirtes, P. (2009) Markov equivalence 
for ancestral graphs. \emph{Annals of Statistics}, 37(5B),2808-2837.
}
\author{
Kayvan Sadeghi
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{MarkEqRcg}}, \code{\link{msep}}
}
\examples{
H1<-matrix(  c(0,100,  0,  0,
	         100,  0,100,  0,
               0,100,  0,100,
               0,  1,100,  0), 4, 4)
H2<-matrix(c(0,0,0,0,1,0,100,0,0,100,0,100,0,1,100,0),4,4)
H3<-matrix(c(0,0,0,0,1,0,0,0,0,1,0,100,0,1,100,0),4,4)
MarkEqMag(H1,H2)
MarkEqMag(H1,H3)
MarkEqMag(H2,H3)
}
\keyword{graphs}
\keyword{Markov equivalence}
\keyword{maximal ancestral graphs}
\keyword{multivariate}