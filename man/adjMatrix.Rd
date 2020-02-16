\name{adjMatrix}
\alias{adjMatrix}
\title{Adjacency matrix of a graph}
\description{
  Transforms the ``edge matrix'' of a graph into
  the adjacency matrix.
}
\usage{
adjMatrix(A)
}
\arguments{
  \item{A}{a square matrix representing the edge matrix of a graph.}
}
\details{
  Given the edge matrix \eqn{A} of a graph, this can be transformed into
  an adjacency matrix \eqn{E} with the formula \eqn{E = (A-I)^T}.
}
\value{
  \item{E}{the adjacency matrix of the graph.}
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{edgematrix}}}
\examples{
amat <- DAG(y ~ x+z, z~u+v)
E <- edgematrix(amat)
adjMatrix(E)
}
\keyword{array}
\keyword{algebra}
\keyword{graphs}
\keyword{multivariate}
