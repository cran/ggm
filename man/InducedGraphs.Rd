\name{InducedGraphs}
\alias{inducedCovGraph}
\alias{inducedConGraph}
\alias{InducedGraphs}
\title{Graphs induced by marginalization or conditioning}
\description{
  Functions to find the induced covariance or concentration graphs after
  conditioning on a set of variables and marginalizing over another set. 
}
\usage{
inducedCovGraph(amat, sel = rownames(amat), cond = NULL)
inducedConGraph(amat, sel = rownames(amat), cond = NULL)
}
\arguments{
  \item{amat}{a square Boolean matrix, the adjacency matrix of a directed
    acyclic graph. The names of rows and of the columns
    are the nodes of the DAG.}
  \item{sel}{a character vector representing a subset
    of selected variables. 
    The elements of the  vector must be a subset of
    the names of the nodes i.e. of \code{rownames(A)}. 
    By default \code{sel} is the set of the nodes of the DAG.
  }
  \item{cond}{a character vector representing the variables on which you
    want to condition. \code{cond} must be disjoint from \code{sel} and
    their union must be a subset of the set of nodes. The set difference
    between the set of nodes and the union of \code{sel} and
    \code{cond} are  the variables over which we marginalize.
    \code{cond} may be the null vector (the default), meaning that you want
    to condition on the empty set.
  }
}
\details{Given a directed acyclic graph representing a set of
  conditional independencies it is possible to obtain other
  graphs of conditional independence implied after marginalizing
  over and conditionig on sets of nodes. Two such graphs are
  the covariance graph and the concentration graph
  (Cox \& Wermuth, 1996, 2003).
}
\value{
  \code{inducedCovGraph} returns the adjacency matrix of the covariance
  graph of the variables in set \code{sel} given the variables
  in set \code{cond}, implied by the original directed acyclic graph
  with adjacency matrix \code{amat}. 
  
  \code{inducedConGraph} returns the adjacency matrix of the concentration
  graph of the variables in set \code{sel} given the variables
  in set \code{cond}, implied by the original directed acyclic graph
  with adjacency matrix \code{amat}.
  
  If \code{sel} is \code{NULL} the functions return 
  the null matrix.
  If \code{cond} is \code{NULL}, the conditioning set is empty and the
  functions return the overall induced covariance or concentration
  matrices of the selected variables.
}
\note{
  If you do not specify \code{sel} you cannot specify a non \code{NULL}
  value of \code{cond}.
  }
  \references{
    Cox, D. R. \& Wermuth, N. (1996). \emph{Multivariate
      dependencies}. London: Chapman \& Hall.
    
    Wermuth, N. \& Cox, D.R. (2003). Joint response graphs
    and separation induced by triangular systems. Submitted and available at
    \url{http://psystat.sowi.uni-mainz.de}.
  }
\author{Giovanni M. Marchetti}
\seealso{\code{\link{DAG}}, \code{\link{UG}},\code{\link{isAcyclic}}}
\examples{
## Define a DAG
dag <- DAG(a ~ x, c ~ b+d, d~ x)
dag
## Induced covariance graph of a, b, d given the empty set.
inducedCovGraph(dag, sel=c("a", "b", "d"), cond=NULL)

## Induced concentration graph of a, b, c given x
inducedConGraph(dag, sel=c("a", "b", "c"), cond="x")

## Overall covariance graph
inducedCovGraph(dag)

## Overall concentration graph
inducedConGraph(dag)

## Induced covariance graph of x, b, d given c, x.
inducedCovGraph(dag, sel=c("a", "b", "d"), cond=c("c", "x"))

## Induced concentration graph of a, x, c given d, b.
inducedConGraph(dag, sel=c("a", "x", "c"), cond=c("d", "b"))

## The DAG on p. 198 of Cox & Wermuth (1996)
dag <- DAG(y1~ y2 + y3, y3 ~ y5, y4 ~ y5)

## Cf. figure 8.7 p. 203 in Cox & Wermuth (1996)
inducedCovGraph(dag, sel=c("y2", "y3", "y4", "y5"), cond="y1")
inducedCovGraph(dag, sel=c("y1", "y2", "y4", "y5"), cond="y3")
inducedCovGraph(dag, sel=c("y1", "y2", "y3", "y4"), cond="y5")

## Cf. figure 8.8 p. 203 in Cox & Wermuth (1996)
inducedConGraph(dag, sel=c("y2", "y3", "y4", "y5"), cond="y1")
inducedConGraph(dag, sel=c("y1", "y2", "y4", "y5"), cond="y3")
inducedConGraph(dag, sel=c("y1", "y2", "y3", "y4"), cond="y5")

## Cf. figure 8.9 p. 204 in Cox & Wermuth (1996)
inducedCovGraph(dag, sel=c("y2", "y3", "y4", "y5"), cond=NULL)
inducedCovGraph(dag, sel=c("y1", "y2", "y4", "y5"), cond=NULL)
inducedCovGraph(dag, sel=c("y1", "y2", "y3", "y4"), cond=NULL)

## Cf. figure 8.10 p. 204 in Cox & Wermuth (1996)
inducedConGraph(dag, sel=c("y2", "y3", "y4", "y5"), cond=NULL)
inducedConGraph(dag, sel=c("y1", "y2", "y4", "y5"), cond=NULL)
inducedConGraph(dag, sel=c("y1", "y2", "y3", "y4"), cond=NULL)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
