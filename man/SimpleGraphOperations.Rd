\name{Simple Graph Operations}
\alias{bd}
\alias{ch}
\alias{pa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simple graph operations}
\description{
Finds the boundary, children, parents of a subset of nodes of a graph.
}
\usage{
bd(nn, A)
ch(nn, A)
pa(nn, A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nn}{a vector of nodes. It may either a numeric vector, or a
    character vector. If it is character vector must be a subset
    of the \code{rownames} of the edge matrix.}
  \item{A}{a square matrix specifying the edge matrix of the graph}
}
\details{For the definitions of these operators see Lauritzen (1996).
}
\value{
  A vector, specifying the boundary or the children or the parents of nodes  \code{nn}
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
bd(2, G)
bd(c(2,3), G)
## find boundary of a subset of nodes of an UG
G <- UG(~ y*x*z + z*h*v)
bd("z", G)
bd(c("y", "x"), G)
bd(5, G)
bd(c(2,5), G)
## children of a subset of nodes of a DAG
G <- DAG(y ~ x+b+a, b~a, x~a)
ch("b", G)
ch(c("b", "x"), G)
ch(2, G)
ch(c(2,3), G)
## parents of a subset of nodes of a DAG
pa("b", G)
pa(c("b", "x"), G)
pa(2, G)
pa(c(2,3), G)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
