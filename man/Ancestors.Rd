\name{Ancestors}
\alias{ancGraph}
\alias{anGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Ancestor graphs}
\description{
Compute the overall ancestors graph from a given DAG.
}
\usage{
ancGraph(A)
anGraph(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{a square Boolean matrix representing
    the edge matrix of the DAG}
}
\details{
  The overall ancestor graph is a graph obtained from a DAG by
  computing the transitive closure. This graph has a non zero
  \eqn{(i,j]} entry iff node \eqn{j} is an ancestor of node \eqn{i}.

  \code{ancGraph} is based on a simple
  formula \eqn{(2 I - A)^{-1}} explained by Wermuth and Cox (2003).

  \code{anGraph} instead uses repeated composition applied to the original
  edge matrix. If \code{A} is not a DAG sometimes \code{ancGraph} fails, while
  \code{anGraph} is correct.
}
\value{
  both functions return the same result, a square Boolean matrix of the
  same order as \code{A} giving the edge matrix of the overall ancestor graph.
  
}
\references{
  Wermuth, N. \& Cox, D.R. (2003). Joint response graphs
  and separation induced by triangular systems. Submitted and available at
  \url{http://psystat.sowi.uni-mainz.de}.
}

\author{Giovanni M. Marchetti}
\seealso{\code{\link{DAG}}, \code{\link{is.acyclic}}}
\examples{
dag <- DAG(y ~ x, x ~ z+u,  u ~ v, a ~ v + b); dag
anGraph(dag) 
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
