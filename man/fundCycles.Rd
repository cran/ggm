\name{fundCycles}
\alias{fundCycles}
\title{Fundamental cycles}
\description{
Finds the list of fundamental cycles of a connected undirected graph.
}
\usage{
fundCycles(amat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{a symmetric matrix with dimnames denoting the adjacency matrix
    of the undirected graph. The graph must be connected, otherwise
    the function returns an error message.}
}
\details{
  All the cycles in an UG can be obtained from combination (ring sum)
  of the set of fundamental cycles.
}
\value{
  a list of matrices with two columns. Every component of the list
  is associated to a cycle. The cycle is described by a
  \eqn{k \times 2} matrix whose rows are the edges of the cycle. If
  there is no cycle the function returns \code{NULL}.
}
\references{
  Thulasiraman, K. & Swamy, M.N.S. (1992). \emph{Graphs: theory and
    algorithms}. New York: Wiley.
}
\author{Giovanni M. Marchetti}
\note{This function is used by \code{cycleMatrix} and \code{isGident}.}
\seealso{\code{\link{UG}},\code{\link{findPath}},
  \code{\link{cycleMatrix}}, \code{\link{isGident}},\code{\link{bfsearch}}}
  \examples{
## Three fundamental cycles
fundCycles(UG(~a*b*d + d*e + e*a*f))
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}

