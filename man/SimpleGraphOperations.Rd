\name{Simple Graph Operations}
\alias{bd}
\alias{ch}
\alias{pa}
\title{Simple graph operations}
\description{
Finds the boundary, children, parents of a subset of nodes of a graph.
}
\usage{
bd(nn, amat)
ch(nn, amat)
pa(nn, amat)
}
\arguments{
  \item{nn}{a vector of nodes. It may either a numeric vector, or a character vector. If it is character vector must be a subset
    of the \code{rownames} of the edge matrix.}
  \item{amat}{a square matrix with dimnames specifying the adjacency matrix of the graph}
}
\details{For definitions of the operators see Lauritzen (1996).
}
\value{
  The operators return a character vector specifying the boundary or the children or the parents of nodes  \code{nn}
  in the graph.  This is a numeric or a character vector depending
  on the mode of \code{nn}.
}
\references{Lauritzen, S. (1996). \emph{Graphical models}. Oxford:
Clarendon Press.}

\author{Giovanni M. Marchetti}
\seealso{\code{\link{UG}}, \code{\link{DAG}}}
\examples{
## find boundary of a subset of nodes of a DAG
G <- DAG(y ~ x+b+a, b~a, x~a)
bd("b", G)
bd(c("b", "x"), G)
bd("x", G)
bd(c("x","b"), G)
## find boundary of a subset of nodes of an UG
G <- UG(~ y*x*z + z*h*v)
bd("z", G)
bd(c("y", "x"), G)
bd("v", G)
bd(c("x","v"), G)
## children of a subset of nodes of a DAG
G <- DAG(y ~ x+b+a, b~a, x~a)
ch("b", G)
ch(c("b", "x"), G)
ch("x", G)
ch(c("a","x"), G)
## parents of a subset of nodes of a DAG
pa("b", G)
pa(c("b", "x"), G)
pa("x", G)
pa(c("x","b"), G)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
