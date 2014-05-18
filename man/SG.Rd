\name{SG}
\alias{SG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summary graph
}
\description{
\code{SG} generates and plots summary graphs after marginalization
 and conditioning.}
\usage{
SG(amat,M=c(),C=c(),showmat=TRUE,plot=FALSE, plotfun = plotGraph, \dots)
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
A subset of the node set of \code{a} that is going to be marginalised over}
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
Sadeghi, K. (2011). Stable classes of graphs containing directed acyclic graphs. 
\emph{Submitted}.

Wermuth, N. (2011). Probability distributions with summary graph structure. 
\emph{Bernoulli}, 17(3),845-879.
}
\author{
Kayvan Sadeghi
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{AG}}, \code{\link{MSG}}, \code{\link{RG}}}
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
M <- c(3,5,6,15,16)
C <- c(4,7)
SG(ex, M, C, plot = TRUE)
SG(ex, M, C, plot = TRUE, plotfun = drawGraph, adjust = FALSE)
}
\keyword{graphs}
\keyword{directed acyclic graph}
\keyword{marginalization and conditioning}
\keyword{summary graph}
