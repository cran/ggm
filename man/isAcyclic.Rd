\name{isAcyclic}
\alias{isAcyclic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graph queries}
\description{
  Checks if a given graph is acyclic.
}
\usage{
isAcyclic(amat)
}
\arguments{
  \item{amat}{a square Boolean matrix with dimnames, the adjacency matrix
    of a graph.}
}
\value{
  a logical value, \code{TRUE} if the graph is acyclic and \code{FALSE}
  otherwise.
}
\references{
  Aho, A.V., Hopcroft, J.E. \& Ullman, J.D. (1983). \emph{Data structures
    and algorithms.} Reading: Addison-Wesley.
}
\author{Giovanni M. Marchetti} 
\examples{
## A cyclic graph
d <- matrix(0,3,3)
rownames(d) <- colnames(d) <- c("x", "y", "z")
d["x","y"] <- d["y", "z"] <- d["z", "x"] <- 1
## Test if the graph is acyclic
isAcyclic(d)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
