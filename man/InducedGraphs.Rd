\name{InducedGraphs}
\alias{inducedCovGraph}
\alias{inducedConGraph}
\alias{InducedGraphs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphs induced by marginalization or conditioning}
\description{
  Functions to find the induced covariance or concentration graphs after
  conditioning on a set of variables and marginalising over another set. 
}
\usage{
inducedCovGraph(A, sel = 1:nrow(A), cond = NULL)
inducedConGraph(A, sel = 1:nrow(A), cond = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{a Boolean matrix with dimnames representing a directed
    acyclic graph. The function does not check if this is the case. See
  the function \code{is.acyclic}}.
  \item{sel}{a vector representing a subset of selected of variables. 
  The vector should be a character vector of the
    names of the variables matching the names of the nodes in
    \code{rownames(A)}. It can be also a numeric vector of indices.
    By default \code{sel} is the set of the nodes of the DAG and
    the \code{cond} is the empty set. 
  }
  \item{cond}{a set of nodes representing the variables on which you
    want to condition. The mode of this vector must match the mode
    of \code{sel}. \code{cond} must be disjoint from \code{sel} and
    their union must be a subset of the set of nodes. The set difference
    between the set of nodes and the union of \code{sel} and
    \code{cond} contains the variables over which we marginalize.
    \code{cond} may be the null vector, meaning that you want
    to condition on the empty set.
  }
}
\details{Given a directed acyclic graph representing a set of
  conditional independencies it is possible to obtain other
  graphs of conditional independence implied after marginalizing
  over and conditionig on sets of nodes. Two such graphs are
  the covariance graph and the concentration graph
  (cfr. Cox \& Wermuth, 1996, 2003).
}
\value{
  \code{inducedCovGraph} returns the edge matrix of the covariance
  graph of the variables in set \code{sel} given the variables
  in set \code{cond}, implied by the original directed acyclic graph
  with edge matrix \code{A}. 
  
  \code{inducedConGraph} returns the edge matrix of the concentration
  graph of the variables in set \code{sel} given the variables
  in set \code{cond}, implied by the original directed acyclic graph
  with edge matrix \code{A}.
  
  If \code{sel} is \code{NULL} the functions return \code{NULL} or the null matrix.
  If \code{cond} is \code{NULL}, the conditioning set is empty and the
  functions return the marginal induced covariance or concentration
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
\seealso{\code{\link{DAG}}, \code{\link{is.acyclic}}}
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

## The variables may be specified by numeric vectors (both)
inducedConGraph(dag, sel=c(1,4,2), cond=c(3, 5))

## The DAG on p. 198 of Cox & Wermuth (1996)
dag <- DAG(y1~ y2 + y3, y3 ~ y5, y4 ~ y5)

## Cfr. figure 8.7 p. 203 in Cox & Wermuth (1996)
inducedCovGraph(dag, sel=c("y2", "y3", "y4", "y5"), cond="y1")
inducedCovGraph(dag, sel=c("y1", "y2", "y4", "y5"), cond="y3")
inducedCovGraph(dag, sel=c("y1", "y2", "y3", "y4"), cond="y5")

## Cfr. figure 8.8 p. 203 in Cox & Wermuth (1996)
inducedConGraph(dag, sel=c("y2", "y3", "y4", "y5"), cond="y1")
inducedConGraph(dag, sel=c("y1", "y2", "y4", "y5"), cond="y3")
inducedConGraph(dag, sel=c("y1", "y2", "y3", "y4"), cond="y5")

## Cfr. figure 8.9 p. 204 in Cox & Wermuth (1996)
inducedCovGraph(dag, sel=c("y2", "y3", "y4", "y5"), cond=NULL)
inducedCovGraph(dag, sel=c("y1", "y2", "y4", "y5"), cond=NULL)
inducedCovGraph(dag, sel=c("y1", "y2", "y3", "y4"), cond=NULL)

## Cfr. figure 8.10 p. 204 in Cox & Wermuth (1996)
inducedConGraph(dag, sel=c("y2", "y3", "y4", "y5"), cond=NULL)
inducedConGraph(dag, sel=c("y1", "y2", "y4", "y5"), cond=NULL)
inducedConGraph(dag, sel=c("y1", "y2", "y3", "y4"), cond=NULL)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
