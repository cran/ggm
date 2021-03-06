\name{makeMG}
\alias{makeMG}
\title{Mixed Graphs}
\description{
  Defines a loopless mixed graph from the directed, undirected and
  undirected components.
}
\usage{
makeMG(dg = NULL, ug = NULL, bg = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dg}{the adjacency matrix of a directed graph specifying
    the arrows of the mixed graph.}
  \item{ug}{the adjacency matrix of an undirected graph specifying the
    lines of the mixed graph.}
  \item{bg}{the adjacency matrix of an undirected graph specifying the
    bidirected edges of the mixed graph.}
}
\details{
  A loopless mixed graph is a mixed graph with three types of edges: undirected,
  directed and bi-directed edges. 
  Note that the three adjacency matrices must have labels and
  may be defined using the functions \code{DG}, \code{DAG} or \code{UG}.  
  The adjacency matrices of the undirected graphs may be just 
  symmetric Boolean matrices.   
}
\value{
  a square matrix obtained by combining the three graph components
  into an adjacency matrix of a mixed graph. The matrix  consists of 
  4 different integers as 
  an \eqn{ij}-element: 0 for a missing 
  edge between \eqn{i} and \eqn{j}, 1 for an arrow from \eqn{i} to \eqn{j}, 
  10 for a full line between 
  \eqn{i} and \eqn{j}, and 100 for a bi-directed arrow between \eqn{i} and \eqn{j}. 
  These numbers are 
  added to be associated with multiple edges of different types. The matrix is 
  symmetric w.r.t full lines and bi-directed arrows.
}
\references{Richardson, T. S. and Spirtes, P. (2002). Ancestral Graph Markov
  Models. \emph{Annals of Statistics}, 30(4), 962--1030.
}
\author{Giovanni M. Marchetti, Mathias Drton}
\seealso{\code{\link{UG}}, \code{\link{DAG}}}
\examples{
## Examples from Richardson and Spirtes (2002)
a1 <- makeMG(dg=DAG(a~b, b~d, d~c), bg=UG(~a*c))  
isAG(a1)    # Not an AG. (a2) p.969    
a2 <- makeMG(dg=DAG(b ~ a, d~c), bg=UG(~a*c+c*b+b*d))           # Fig. 3 (b1) p.969  
isAG(a1)
a3 <- makeMG(ug = UG(~ a*c), dg=DAG(b ~ a, d~c), bg=UG(~ b*d)) # Fig. 3 (b2) p.969
a5 <- makeMG(bg=UG(~alpha*beta+gamma*delta), dg=DAG(alpha~gamma,
delta~beta))  # Fig. 6 p. 973
## Another Example
a4 <- makeMG(ug=UG(~y0*y1), dg=DAG(y4~y2, y2~y1), bg=UG(~y2*y3+y3*y4))  
## A mixed graphs with double edges. 
mg <- makeMG(dg = DG(Y ~ X, Z~W, W~Z, Q~X), ug = UG(~X*Q), 
bg = UG(~ Y*X+X*Q+Q*W + Y*Z) )
## Chronic pain data: a regression graph
chronic.pain <- makeMG(dg = DAG(Y ~ Za, Za ~ Zb + A, Xa ~ Xb, 
Xb ~ U+V, U ~ A + V, Zb ~ B, A ~ B), bg = UG(~Za*Xa + Zb*Xb))
}
\keyword{graphs} 
\keyword{ancestral graph}      
\keyword{mixed graph}
\keyword{models}
\keyword{multivariate}
