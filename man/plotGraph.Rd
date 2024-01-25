\name{plotGraph}
\alias{plotGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot of a mixed graph}
\description{
Plots a mixed graph from an adjacency matrix,
a \code{graphNEL} object, an \code{\link{igraph}} object, or
a descriptive vector.
}
\usage{
plotGraph(a, dashed = FALSE, tcltk = TRUE, layout = layout.auto,
directed = FALSE, noframe = FALSE, nodesize = 15, vld = 0, vc = "gray",
vfc = "black", colbid = "FireBrick3", coloth = "black", cex = 1.5, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
An adjacency matrix: a matrix that consists of 4 different integers as an \eqn{ij}-element:
0 for a missing edge between \eqn{i} and \eqn{j}, 1 for an arrow from \eqn{i} to \eqn{j},
10 for a full line between
\eqn{i} and \eqn{j}, and 100 for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers can be added to generate multiple edges of different types. The matrix must be symmetric w.r.t full lines and bi-directed arrows.
Or a graph that can be a \code{graphNEL} or an \code{\link{igraph}} object.Or a vector of length \eqn{3e}, where \eqn{e} is the number of edges of the graph, that is a sequence of triples (type,node1label,node2label). The type of edge can be \code{"a"} (arrows from node1 to node2), \code{"b"} (arcs), and \code{"l"} (lines).
}
\item{dashed}{A logical value. If \code{TRUE} the bi-directed edges are plotted as undirected dashed edges.}
\item{tcltk}{A logical value. If \code{TRUE} the function opens a tcltk device to plot the graphs, allowing the interactive manimulation of the graph. If \code{FALSE}the function opens a standard device without interaction.}
\item{layout}{The name of a function used to compute the (initial) layout of the graph. The default is \code{layout.auto}. This can be further adjusted if \code{tcltk} is \code{TRUE}.}
\item{directed}{A logical value. If \code{FALSE} a symmetric adjacency matrix with entries 1 is interpreted as an undirected graph. If \code{TRUE} it is interpreted as a directed graph with double arrows. If \code{a} is not an adjacency matrix, it is ignored.}
\item{noframe}{A logical value. If \code{TRUE}, then the nodes are not circled.}
\item{nodesize}{An integer denoting the size of the nodes (default 15). It can be
increased to accommodate larger labels.}
\item{vld}{An integer defining the distance between a vertex and its label. Defaults to 0.}
\item{vc}{Vertex color. Default is "gray".}
\item{vfc}{Vertex frame color. Default is "black".}
\item{colbid}{Color of the bi-directed edges. Default is "FireBrick3".}
\item{coloth}{Color of all the other edges. Default is "black".}
\item{cex}{An integer (defaults to 1) to adjust the scaling of the font of the labels.}
\item{\dots}{Further arguments to be passed to \code{plot} or \code{tkplot}.}
}
\details{
\code{plotGraph} uses \code{\link{plot}} and \code{\link{tkplot}} in \pkg{\link{igraph}} package.
}
\value{
Plot of the associated graph and returns invisibly a list with two slots:
\code{tkp.id},
\code{graph}, the input graph as an \code{igraph} object.
The id can be used to get the layout of the adjusted graph.
The bi-directed edges are plotted in red.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Kayvan Sadeghi, Giovanni M. Marchetti
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{grMAT}}, \code{\link{tkplot}}, \code{\link{drawGraph}},
\code{\link{plot.igraph}}
}
\examples{
exvec<-c("b",1,2,"b",1,14,"a",9,8,"l",9,11,
         "a",10,8,"a",11,2,"a",11,9,"a",11,10,
         "a",12,1,"b",12,14,"a",13,10,"a",13,12)
plotGraph(exvec)
############################################
amat<-matrix(c(0,11,0,0,10,0,100,0,0,100,0,1,0,0,1,0),4,4)
plotGraph(amat)
plotGraph(makeMG(bg = UG(~a*b*c+ c*d), dg = DAG(a ~ x + z, b ~ z )))
plotGraph(makeMG(bg = UG(~a*b*c+ c*d), dg = DAG(a ~ x + z, b ~ z )), dashed = TRUE)
# A graph with double and triple edges
G <-
structure(c(0, 101, 0, 0, 100, 0, 100, 100, 0, 100, 0, 100, 0,
111, 100, 0), .Dim = c(4L, 4L), .Dimnames = list(c("X", "Z",
"Y", "W"), c("X", "Z", "Y", "W")))
plotGraph(G)
# A regression chain graph with longer labels
 plotGraph(makeMG(bg = UG(~Love*Constraints+ Constraints*Reversal+ Abuse*Distress),
   dg = DAG(Love ~ Abuse + Distress, Constraints ~ Distress, Reversal ~ Distress,
   Abuse ~ Fstatus, Distress ~ Fstatus),
   ug = UG(~Fstatus*Schooling+ Schooling*Age)),
   dashed = TRUE, noframe = TRUE)
# A graph with 4 edges between two nodes.
G4 = matrix(0, 2, 2); G4[1,2] = 111; G4[2,1] = 111
plotGraph(G4)
}
\keyword{graphs}
\keyword{adjacency matrix}
\keyword{mixed graphs}
\keyword{plot}
