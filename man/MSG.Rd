\name{MSG}
\alias{MSG}
\title{
Maximal summary graph
}
\description{
\code{MAG} generates and plots maximal summary graphs after marginalization
 and conditioning.}
\usage{
MSG(amat,M=c(),C=c(),showmat=TRUE,plot=FALSE, plotfun = plotGraph, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{amat}{
An adjacency matrix of a MAG, or a graph that can be a \code{graphNEL} or an \code{\link{igraph}} object
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
%  \item{byrownames}{
%A logical value. \code{TRUE} (by default) if \code{M} and \code{C} are
%the subsets of the labels of the node set of the
%graph or the row-name of the matrix of the first argument. \code{FALSE}
%if \code{M} and \code{C} are the subsets of the set \eqn{\{1, \dots, n\}}, where \eqn{n} is the
%number of the nodes.
%}
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
This function uses the functions \code{\link{SG}} and \code{\link{Max}}.
}
\value{
A matrix that consists 4 different integers as an \eqn{ij}-element: 0 for a missing
edge between \eqn{i} and \eqn{j}, 1 for an arrow from \eqn{i} to \eqn{j}, 10 for a full line between
\eqn{i} and \eqn{j}, and 100 for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers are
added to be associated with multiple edges of different types. The matrix is
symmetric w.r.t full lines and bi-directed arrows.
}
\references{
Richardson, T.S. and Spirtes, P. (2002). Ancestral graph Markov models. \emph{Annals
of Statistics}, 30(4), 962-1030.


Sadeghi, K. (2013). Stable mixed graphs.
\emph{Bernoulli} 19(5B), 2330–2358.

Sadeghi, K. and Lauritzen, S.L. (2014). Markov properties for loopless mixed graphs. \emph{Bernoulli} 20(2), 676-696.

Wermuth, N. (2011). Probability distributions with summary graph structure.
\emph{Bernoulli}, 17(3), 845-879.}
\author{
Kayvan Sadeghi
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{MAG}}, \code{\link{Max}}, \code{\link{MRG}}, \code{\link{SG}}}
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
             0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0), 16, 16, byrow=TRUE)
M <- c(3,5,6,15,16)
C <- c(4,7)
MSG(ex,M,C,plot=TRUE)
###################################################
H<-matrix(c(0,100,1,0,100,0,100,0,0,100,0,100,0,1,100,0),4,4)
Max(H)
}
\keyword{graphs}
\keyword{directed acyclic graph}
\keyword{marginalisation and conditioning}
\keyword{maximality of graphs}
\keyword{summary graph}
