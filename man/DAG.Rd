\name{DAG}
\alias{DAG}
\title{Directed acyclic graphs (DAGs)}
\description{
  A simple way to define a DAG by means of regression model
  formulae.
  }
\usage{
DAG(\dots, order = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{a sequence of model formulae}
  \item{order}{logical, defaulting to \code{FALSE}. If \code{TRUE} the nodes of the
    DAG are permuted according to the topological order. If
    \code{FALSE} the nodes are in the order they first appear in the model
    formulae (from left to right).}
}
\details{
  The DAG is defined by a sequence of recursive regression models.
  Each regression is defined by a model formula.
  For each formula the response defines a node of the graph and
  the explanatory variables the parents of that node. If the
  regressions are not recursive the function returns an error message.

  Some authors prefer the terminology acyclic directed graphs (ADG).
}
\value{the adjacency matrix of the DAG, i.e. 
  a square Boolean matrix of order equal to the number of nodes of the
  graph and a one in position \eqn{(i,j)} if there is an arrow from
  \eqn{i} to \eqn{j} and zero otherwise.  The rownames of the adjacency
  matrix are the nodes of the DAG. 
  
  If \code{order = TRUE} the
  adjacency matrix is permuted to have  parents before  children.
  This can always be done (in more than one way) for DAGs. The resulting
  adjacency matrix is upper triangular. 
}
\references{
  Lauritzen, S. (1996). \emph{Graphical models}. Oxford:
  Clarendon Press.
}
\author{G. M. Marchetti}
\note{The model formulae may contain interactions, but they are
  ignored in the graph.} 
\seealso{\code{\link{UG}}, \code{\link{topSort}}, \code{\link{edgematrix}}, \code{\link{fitDag}}}
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

## The order can be changed
DAG(y ~ z, y ~ x + u + v,  u ~ v + z)

## If you want to order the nodes (topological sort of the DAG)
DAG(y ~ z, y ~ x + u + v,  u ~ v + z, order=TRUE)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
