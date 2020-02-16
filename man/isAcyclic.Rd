\name{isAcyclic}
\alias{isAcyclic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graph queries}
\description{
  Checks if a given graph is acyclic.
}
\usage{
isAcyclic(amat, method = 2)
}
\arguments{
  \item{amat}{a square Boolean matrix with dimnames, the adjacency matrix of a graph.}
  \item{method}{an integer 1 or 2 specifying the method used. If \code{method=1} the  function calls the function \code{clusters} in package \code{igraph} to find the strong components: two nodes v and w are in the same strong component iff there are directed paths from v to w and from w to v. If \code{method=2} the function uses the \code{ggm} function \code{transClos}. Method 1 is faster.}
}
\value{
  a logical value, \code{TRUE} if the graph is acyclic and \code{FALSE}
  otherwise.
}
\references{
  Aho, A.V., Hopcroft, J.E. \& Ullman, J.D. (1983). \emph{Data structures
    and algorithms.} Reading: Addison-Wesley.
}
\author{David Edwards, Giovanni M. Marchetti} 
\examples{
## A cyclic graph
d <- matrix(0,3,3)
rownames(d) <- colnames(d) <- c("x", "y", "z")
d["x","y"] <- d["y", "z"] <- d["z", "x"] <- 1
## Test if the graph is acyclic
isAcyclic(d)
isAcyclic(d, method = 1)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
