\name{DAG}
\alias{DAG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Defining directed acyclic graphs (DAGs)}
\description{
  A simple way to define a DAG by means of regression model
  formulae.
  }
\usage{
DAG(..., order = FALSE, test = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{a sequence of model formulae}
  \item{order}{logical, defaulting to \code{FALSE}. If \code{TRUE} the nodes of the
    DAG are permuted to obtain and upper triangular edge matrix. If
    \code{FALSE} the nodes are in the order they first appear in the model
    formulae (from left to right).}
  \item{test}{logical value. If \code{TRUE} (the default) then it is
    tested if the graph is acyclic. Else if \code{test=FALSE} the
    test is bypassed. In this case the function may return a cyclic
    graph.}  
}
\details{
  The DAG is defined by a sequence of recursive regression models.
  Each regression is defined by a model formula.
  For each formula the response defines a node of the graph and
  the explanatory variables the parents of that node. If the
  regressions are not recursive the function returns an error message
  (if \code{test=TRUE}).
  
  The function returns the edge matrix of the graph, i.e. 
  a square Boolean matrix of order equal to the number of nodes of the
  graph and a one in position \eqn{(i,j)} if there is an arrow from
  \eqn{j} to \eqn{i} and zero otherwise. By default this matrix
  has ones along the main diagonal. Optionally (if \code{order = TRUE}) the
  edge matrix is permuted to be in upper triangular form (this can
  always be done for DAGs). The rownames of the edge matrix are
  the nodes of the DAG.
}
\value{a square Boolean matrix, with dimnames, the edge matrix of the DAG.}
\references{
  Lauritzen, S. (1996). \emph{Graphical models}. Oxford:
  Clarendon Press.
  
  Wermuth, N. \& Cox, D.R. (2003). Joint response graphs
  and separation induced by triangular systems. Submitted and available at
  \url{http://psystat.sowi.uni-mainz.de}.
}
\author{G. M. Marchetti}
\note{The model formulae may contain interactions, but they are
  ignored in the graph.} 
\seealso{\code{\link{fitDag}}, \code{\link{UG}}, \code{\link{topSort}}}
\examples{
## A Markov chain
DAG(y ~ x, x ~ z, z ~ u)

## Another DAG
DAG(y ~ x + z + u, x ~ u, z ~ u)

## A DAG with an isolated node
DAG(v ~ v, y ~ x + z, z ~ w + u)

## There can be repetitions
DAG(y ~ x + u + v, y ~ z, u ~ v + z)

## Interactions are ignored
DAG(y ~ x*z + z*v, x ~ z)

## A cyclic graph returns an error!
\dontrun{DAG(y ~ x, x ~ z, z ~ y)}

## Forcing a cyclic graph!
DAG(y ~ x, x ~ z, z ~ y, test=FALSE)

## The order can be changed
DAG(y ~ z, y ~ x + u + v,  u ~ v + z)

## If you want to order the nodes (top sort of the DAG)
DAG(y ~ z, y ~ x + u + v,  u ~ v + z, order=TRUE)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
