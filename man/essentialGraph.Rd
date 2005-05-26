\name{essentialGraph}
\alias{essentialGraph}
\title{Essential graph}
\description{
  Find the essential graph from a given directed acyclic graph. 
}
\usage{
essentialGraph(dagx)
}
\arguments{
  \item{dagx}{a square binary matrix, the adjacency matrix of a directed
    acyclic graph. The names of rows and of the columns
    are the nodes of the DAG.}
}
\details{Converts a DAG into the Essential Graph.
  Is implemented by the algorithm by D.M.Chickering (1995).
}
\value{
  returns the adjacency matrix of the essential graph.
}
\references{
  Chickering, D.M. (1995). A transformational characterization of
  equivalent Bayesian network structures. \emph{Proceedings of Eleventh
    Conference on Uncertainty in Artificial Intelligence}, Montreal, QU,
  87-98. Morgan Kaufmann.
  
  \url{http://research.microsoft.com/~dmax/publications/uai95.pdf} 
  }
  \author{Giovanni M. Marchetti, translating a MATLAB function
    by Tomas Kocka, AAU}
\seealso{\code{\link{DAG}}, \code{\link{InducedGraphs}}}
\examples{
dag = DAG(U ~ Y+Z, Y~X, Z~X)
essentialGraph(dag)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
