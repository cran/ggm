\name{edgematrix}
\alias{edgematrix}
\title{Edge matrix of a graph}
\description{
  Transforms the adjacency matrix of a graph into
  an ``edge matrix''.
}
\usage{
edgematrix(E, inv=FALSE)
}
\arguments{
  \item{E}{a square matrix, representing the adjacency matrix of a
    graph.}
  \item{inv}{a logical value.}
}
\details{
  In some matrix computations for graph objects the adjacency matrix
  of the graph is transformed into an ``edge matrix''. Briefly,
  if \eqn{E} is the adjacency matrix of the
  graph, the edge matrix is \eqn{A = sign(E+I)^T=[a_{ij}]}. 
  Thus, \eqn{A} has ones along the diagonal
  and if the graph has no edge between nodes \eqn{i} and \eqn{j} the entries
  \eqn{a_{i,j}} and \eqn{a_{j,i}} are both zero.  
  If there is an arrow from \eqn{j} to \eqn{i} 
  \eqn{a_{i,j}=1} and  \eqn{a_{j,i} = 0}. If there is an undirected edge, both
  \eqn{a_{i,j}=a_{j,i}=1}.  
}
\value{
  \item{A}{the edge matrix of the graph.
   If \code{TRUE} the nodes are sorted in
  inverted topological order and the edge matrix is upper triangular.}
}
\references{Wermuth, N. (2003). Analysing social science data with
  graphical Markov models. In: \emph{Highly Structured Stochastic
    Systems.} P. Green, N. Hjort \& T. Richardson (eds.),
  47--52. Oxford: Oxford University Press.}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{adjMatrix}}}
\examples{
amat <- DAG(y ~ x+z, z~u+v)
amat
edgematrix(amat)
edgematrix(amat, inv=TRUE)
}
\keyword{array}
\keyword{algebra}
\keyword{graphs}
\keyword{multivariate}
