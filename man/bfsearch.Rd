\name{bfsearch}
\alias{bfsearch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Breadth first search}
\description{Breadth-first search of a connected undirected graph.
}
\usage{
bfsearch(amat, v = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{a symmetric matrix with dimnames specifying the adjacency
    matrix of the undirected graph}
  \item{v}{an integer, indicating the starting node of the
    search. Defaults to the first node.
  }
}
\details{
  Breadth-first search is a systematic method for exploring a graph.
  The algorithm is taken from Aho, Hopcroft \& Ullman (1983).
}
\value{
  \item{tree}{the edge matrix of the resulting spanning tree}
  \item{branches}{a matrix with two columns, giving the indices of the
    branches of the spanning tree}
  \item{chords}{a matrix with two columns, giving the indices of the
    chords of the spanning tree}
}
\references{
Aho, A.V., Hopcrtoft, J.E. \& Ullman, J.D. (1983). \emph{Data structures
  and algorithms.} Reading: Addison-Wesley.

  Thulasiraman, K. \& Swamy, M.N.S. (1992). \emph{Graphs: theory and
    algorithms}. New York: Wiley.
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{UG}}, \code{\link{findPath}}, \code{\link{cycleMatrix}}}
\examples{
## Finding a spanning tree of the butterfly graph
bfsearch(UG(~ a*b*o + o*u*j))
## Starting from another node
bfsearch(UG(~ a*b*o + o*u*j), v=3)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
