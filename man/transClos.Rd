\name{transClos}
\alias{transClos}
\title{Transitive closure of a graph}
\description{
  Computes the transitive closure of a graph (undirected or directed
  acyclic).}
\usage{
transClos(amat)
}
\arguments{
  \item{amat}{a Boolean matrix with dimnames representing the adjacency
    matrix
    of a graph.}
}
\details{
The transitive closure of a directed graph with adjacency matrix \eqn{A} is a
graph with adjacency matrix \eqn{A^*} such that \eqn{A^*_{i,j} = 1}
if there is a directed path from \eqn{i} to \eqn{j}. The transitive
closure of an undirected graph is defined similarly (by substituting
path to directed path).
}
\value{
  \item{A}{The adjacency matrix of the transitive closure.}
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{DAG}}, \code{\link{UG}}}
\examples{
## Closure of a DAG
d <- DAG(y ~ x, x ~ z)
transClos(d)

## Closure of an UG
g <- UG(~ x*y*z+z*u+u*v)
transClos(g)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
