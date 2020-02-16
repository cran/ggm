\name{cmpGraph}
\alias{cmpGraph}
\title{The complementary graph}
\description{
Finds the complementary graph of an undirected graph.
}
\usage{
cmpGraph(amat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{the adjacency matrix of an undirected graph}
}
\details{
  The complementary graph of an UG is the graph that has the same set
  of nodes and an undirected edge connecting \eqn{i} and \eqn{j}
  whenever there is not an \eqn{(i,j)} edge in the original UG.
}
\value{
the edge matrix of the complementary graph.
}
\references{Lauritzen, S. (1996). \emph{Graphical models}. Oxford:
Clarendon Press.}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{UG}}, \code{\link{DAG}}}
\examples{
## A chordless four-cycle
four <- UG(~ a*b + b*d + d*e + e*a)
four
cmpGraph(four)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
