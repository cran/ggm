\name{MAG}
\alias{MAG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximal ancestral graph 
}
\description{
\code{MAG} generates and plots maximal ancestral graphs after marginalisation
 and conditioning.}
\usage{
MAG(amat,M=c(),C=c(),showmat=TRUE,plot=FALSE, plotfun = plotGraph, \dots)   
}
\arguments{
\item{amat}{
An adjacency matrix, or a graph that can be a \code{graphNEL} or an \code{\link{igraph}} object 
or a vector of length \eqn{3e}, where \eqn{e} is the number of edges of the graph, 
that is a sequence of triples (type, node1label, node2label). The type 
of edge can be \code{"a"} (arrows from node1 to node2), \code{"b"} (arcs), and 
\code{"l"} (lines).
}
  \item{M}{
A subset of the node set of \code{a} that is going to be marginalized over}
  \item{C}{
Another disjoint subset of the node set of \code{a} that is going to be
conditioned on.
}
%\item{byrownames}{
%A logical value. \code{TRUE} (by default) if \code{M} and \code{C} are 
%the subsets of the labels of the node set of the 
%graph or the row-name of the matrix of the first argument. \code{FALSE} 
%if \code{M} and \code{C} are the subsets of the set \eqn{\{1,\dots,n\}}, where \eqn{n} is the 
%number of the nodes.}
\item{showmat}{
A logical value. \code{TRUE} (by default) to print the generated matrix.
}
\item{plot}{
A logical value, \code{FALSE} (by default). \code{TRUE} to plot 
the generated graph. 
}
\item{plotfun}{
Function to plot the graph when \code{plot == TRUE}. Can be \code{plotGraph} (the default) or \code{drawGraph}.
}
\item{\dots}{Further arguments passed to \code{plotfun}.}
}
\details{
This function uses the functions \code{\link{AG}} and \code{\link{Max}}.
}
\value{
A matrix that consists 4 different integers as an \eqn{ij}-element: 0 for a missing 
edge between \eqn{i} and \eqn{j}, 1 for an arrow from \eqn{i} to \eqn{j}, 10 for a full line between 
\eqn{i} and \eqn{j}, and 100 for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers are 
added to be associated with multiple edges of different types. The matrix is 
symmetric w.r.t full lines and bi-directed arrows.
}
\references{
Richardson, T. S. and Spirtes, P. (2002). Ancestral graph Markov models. \emph{Annals 
of Statistics}, 30(4), 962-1030.

Sadeghi, K. (2011). Stable classes of graphs containing directed acyclic graphs. 
\emph{Submitted}.

Sadeghi, K.  and Lauritzen, S.L. (2011). Markov properties for loopless mixed graphs.
\emph{Submitted}. URL \url{http://arxiv.org/abs/1109.5909}.}
\author{
Kayvan Sadeghi
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{AG}}, \code{\link{Max}}, \code{\link{MRG}}, \code{\link{MSG}}
}
\examples{
ex<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ##The adjacency matrix of a DAG
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
             0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0), 16, 16, byrow = TRUE)
M <- c(3,5,6,15,16)
C <- c(4,7)
MAG(ex, M, C, plot=TRUE)
###################################################
H <- matrix(c(0,100,1,0,100,0,100,0,0,100,0,100,0,1,100,0),4,4)
Max(H)
}

\keyword{ancestral graph}
\keyword{directed acyclic graph}
\keyword{marginalization and conditioning}
\keyword{maximality of graphs}
