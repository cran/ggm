\name{topSort}
\alias{topSort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Topological sort}
\description{
Permutates the nodes of a directed acyclic graph in such a way its edge
matrix is upper triangular.
}
\usage{
topSort(gmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gmat}{a square Boolean matrix with dimnames representing the
    edge matrix of a directed acyclic graph.}
}
\details{
The edge matrix of a directed acyclic graph (DAG) can always transformed
into an upper triangular matrix by a permutation applied both to rows
and columns.  The ordering needs not to be unique.
The algorithm to perform the ordering of the nodes is
called topological sort.

}
\value{
  \item{triu}{the sorted edge matrix in upper triangular form. If the
    edge matrix is already in upper triangular form it is leaved untrasformed.}
  \item{perm}{the vector of the permutation.}
}
\note{
  The order of the nodes defined by \code{DAG} is that of their first
  appearance in the model formulae (from left to right).
}
\references{
  Aho, A.V., Hopcrtoft, J.E. \& Ullman, J.D. (1983). \emph{Data structures
    and algorithms.} Reading: Addison-Wesley.
    
    Lauritzen, S. (1996). \emph{Graphical models}. Oxford:
Clarendon Press.
}
\author{Giovanni M. Marchetti} 
\seealso{\code{\link{DAG}}, \code{\link{is.acyclic}}}
\examples{
## A simple example
dag <- DAG(a ~ b, c ~ a + b, d ~ c + b)
dag
m <- topSort(dag)
m
dag[m$perm, m$perm]
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
