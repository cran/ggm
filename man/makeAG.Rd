\name{makeAG}
\alias{makeAG}
\title{Ancestral Graphs}
\description{
  Defines an ancestral graph from the directed, undirected and
  undirected components and checks if the components are
  compatible with an ancestral graph.
}
\usage{
makeAG(dag = NULL, ug = NULL, bg = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dag}{the adjacency matrix of a directed acyclic graph specifying
    the arrows of the ancestral graph.}
  \item{ug}{the adjacency matrix of an undirected graph specifying the
    lines of the ancestral graph.}
  \item{bg}{the adjacency matrix of an undirected graph specifying the
    bidirected edges of the ancestral graph.}
}
\details{
  An ancestral graph is a mixed graph with three types of edges: undirected,
  directed and bidirected edges. The following conditions must hold:
  (i) no undirected edge meets an arrowhead; (ii) no directed cycles;
  (iii) spouses cannot be ancestors. For details see Richardson \&
  Spirtes (2002).

  The function checks if, given the matrices
  of the undirected, directed and bidirected edges, the above three conditions
  are respected. If so, a resulting adjacency matrix \eqn{E=[e_{ij}} is
  returned, with the following convention. If \eqn{(i,j)} is a directed
  edge, then \eqn{e_{ij} = 1} and \eqn{e_{ji}=0}.
  If \eqn{(i,j)} is an undirected
  edge, then \eqn{e_{ij} = e_{ji} = 1}. Finally, if \eqn{(i,j)} is a bidirected
  edge, then \eqn{e_{ij} = e_{ji} = 2}.

  Note that the three adjacency matrices must have labels and
  may be defined using the functions \code{DAG} and \code{UG}.
}
\value{
  a square matrix obtained by combining the three graph components
  into an adjacency matrix of an ancestral graph.
  See the details for the coding of the adjacency matrix.
}
\references{Richardson, T. Spirtes, P. (2002). Ancestral Graph Markov
  Models. \emph{Annals of Statistics}. 30, 4, 962--1030.
}
\author{Giovanni M. Marchetti, Mathias Drton}
\seealso{\code{\link{UG}}, \code{\link{DAG}}}
\examples{
## Examples from Richardson and Spirtes (2002)
\dontrun{a1 <- makeAG(dag=DAG(a~b, b~d, d~c), bg=UG(~a*c))}  # Not an AG. (a2) p.969
a2 <- makeAG(dag=DAG(b ~ a, d~c), bg=UG(~a*c+c*b+b*d))           # Fig. 3 (b1) p.969
 a3 <- makeAG(ug = UG(~ a*c), dag=DAG(b ~ a, d~c), bg=UG(~ b*d)) # Fig. 3 (b2) p.969
 a5 <- makeAG(bg=UG(~alpha*beta+gamma*delta), dag=DAG(alpha~gamma,
delta~beta))  # Fig. 6 p. 973
## Another Example
 a4 <- makeAG(ug=UG(~y0*y1), dag=DAG(y4~y2, y2~y1), bg=UG(~y2*y3+y3*y4))
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
