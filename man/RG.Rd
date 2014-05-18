\name{RG}
\alias{RG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ribbonless graph
}
\description{
\code{RG} generates and plots ribbonless graphs (a modification of MC graph
to use m-separation) after marginalization and conditioning.}
\usage{
RG(amat,M=c(),C=c(),showmat=TRUE,plot=FALSE, plotfun = plotGraph, \dots)
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
\value{
A matrix that consists 4 different integers as an \eqn{ij}-element: 0 for a missing 
edge between \eqn{i} and \eqn{j}, 1 for an arrow from \eqn{i} to \eqn{j}, 10 for a full line between 
\eqn{i} and \eqn{j}, and 100 for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers are 
added to be associated with multiple edges of different types. The matrix is 
symmetric w.r.t full lines and bi-directed arrows.
}
\references{
Koster, J.T.A. (2002). Marginalizing and conditioning in graphical models. 
\emph{Bernoulli}, 8(6), 817-840.

Sadeghi, K. (2011). Stable classes of graphs containing directed acyclic graphs. 
\emph{Submitted}.
}
\author{
Kayvan Sadeghi
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{AG}},, \code{\link{MRG}}, \code{\link{SG}}
}
\examples{
	ex <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ##The adjacency matrix of a DAG
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
	               0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0),16,16, byrow = TRUE)

M<-c(3,5,6,15,16)
C<-c(4,7)
RG(ex,M,C,plot=TRUE)
}
\keyword{graphs}
\keyword{directed acyclic graph}
\keyword{marginalisation and conditioning}
\keyword{MC graph}
\keyword{ribbonless graph}
