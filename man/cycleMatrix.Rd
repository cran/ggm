\name{cycleMatrix}
\alias{cycleMatrix}
\title{Fundamental cycles}
\description{
Finds the matrix of fundamental cycles of a connected undirected graph.
}
\usage{
cycleMatrix(amat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{a symmetric matrix with dimnames denoting the adjacency matrix
    of the undirected graph. The graph must be connected, otherwise
    the function returns an error message.}
}
\details{
  All the cycles in an UG can be obtained from combination (ring sum)
  of the set of fundamental cycles. The matrix of fundamental cycles
  is a Boolean matrix having as rows the fundamental cycles and as
  columns the edges of the graph. If an entry is one then the edge
  associated to that column belongs to the cycle associated to the row.
}
\value{
  a Boolean matrix of the fundamental cycles of the undirected graph.
  If   there is no cycle the function returns \code{NULL}.}
\references{
  Thulasiraman, K. & Swamy, M.N.S. (1992). \emph{Graphs: theory and
    algorithms}. New York: Wiley.
}
\author{Giovanni M. Marchetti}
\note{This function is used by \code{isGident}. The row sum of the
  matrix gives the length of the cycles.}
\seealso{\code{\link{UG}}, \code{\link{findPath}},
  \code{\link{fundCycles}}, \code{\link{isGident}}, \code{\link{bfsearch}}}
  \examples{
## Three cycles
cycleMatrix(UG(~a*b*d+d*e+e*a*f))
## No cycle
 cycleMatrix(UG(~a*b))
## two cycles: the first is even and the second is odd
cm <- cycleMatrix(UG(~a*b+b*c+c*d+d*a+a*u*v))
apply(cm, 1, sum)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
