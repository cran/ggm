\name{unmakeMG}
\alias{unmakeMG}
\title{Loopless mixed graphs components}
\description{
  Splits the adjacency matrix of a loopless mixed graph into  three components:
  directed, undirected and bi-directed.
}
\usage{
unmakeMG(amat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{a square matrix, with dimnames, representing a 
  loopless mixed graph. The matrix  consists of 4 different integers as 
an \eqn{ij}-element: 0 for a missing 
edge between \eqn{i} and \eqn{j}, 1 for an arrow from \eqn{i} to \eqn{j}, 10 for a full line between 
\eqn{i} and \eqn{j}, and 100 for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers are 
added to be associated with multiple edges of different types. The matrix is 
symmetric w.r.t full lines and bi-directed arrows.
}
}
\value{
  It is the inverse of \code{makeAG}. It returns the following
  components.
  \item{dg}{the adjacency matrix of the directed edges.}
  \item{ug}{the adjacency matrix of the undirected edges.}
  \item{bg}{the adjacency matrix of the bi-directed edges.}
} 
\details{The matrices \code{ug}, and \code{bg} are just symmetric Boolean
matrices.}
\author{Mathias Drton, Giovanni M. Marchetti}
\seealso{\code{\link{makeMG}}}
\examples{
ag <- makeMG(ug=UG(~y0*y1), dg=DAG(y4~y2, y2~y1), bg=UG(~y2*y3+y3*y4))  
isAG(ag)
unmakeMG(ag)
}
\keyword{graphs}    
\keyword{ancestral graph}      
\keyword{mixed graph}
\keyword{models}
\keyword{multivariate}
