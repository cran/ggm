\name{drawGraph}
\alias{drawGraph}
\title{Drawing a graph with a simple point and click interface.}
\description{Draw a graph from its adjacency matrix representation.}
\usage{
drawGraph(amat, coor = NULL, adjust = FALSE, alpha = 1.5,
                beta = 3, lwd = 1, ecol = "blue", bda = 0.1, layout = layout.auto)
}
\arguments{
  \item{amat}{the adjacency matrix representation of the graph. This
    can be an undirected graph, a directed acyclic graph or a mixed 
    graph with at most a summary graph structure. See also \code{\link{plotGraph}}}.
  \item{coor}{an optional matrix of dimensions \code{p} x 2 where
    \eqn{p} is the number of vertices of the graph. If \code{coor=NULL}
    then the function chooses a default position for the nodes.}
  \item{adjust}{a logical value, defaults to \code{FALSE}. If \code{TRUE} the graph is plotted and the system waits until the mouse button is pressed (same behaviour
    of \code{locator} function.}  
  \item{alpha}{a positive value between controlling the
    distance from the end of the edges to the nodes of the graph.}
  \item{beta}{a positive value controlling the distance of the labels of
    the variables from the nodes.}
  \item{lwd}{line width of the edges (default: 1).}
  \item{ecol}{color of the edges (default: "blue").}
  \item{bda}{bidirected edge arrow length (default: 0.1).}
  \item{layout}{The name of a function used to compute the (initial) layout of the graph. The default is \code{layout.auto}. This can be further adjusted if \code{adjust} is \code{TRUE}.}
}
\details{
  The function is a very simple tool useful for displaying small graphs,
  with a rudimentary interface for moving nodes and edges of a given
  graph and  adjusting the final plot. For better displays use
  \pkg{dynamicGraph} or \pkg{Rgraphviz} package in Bioconductor project.
}
\value{
  The function plots the graph with a initial positioning of the nodes,
  as  specified by  \code{coor} and remains in a waiting state.
  The position of each node can be shifted
  by pointing and clicking (with the first mouse button) close to the node.
  When the mouse button is pressed the node which is closer to the
  selected point is moved to that position.
  Thus, one must be careful to click closer to
  the selected node than to any other node. 
  The nodes can be moved to any position by repeating the
  previous operation.
  The adjustment process is terminated by pressing any mouse
  button other than the first.

  At the end of the process, the function returns invisibly the
  coordinates of the nodes. The coordinates may be used later to
  redisplay the graph.
}
\references{
\pkg{dynamicGraph}, \pkg{Rgraphwiz}, \url{http://www.bioconductor.org}.
  
  GraphViz, Graph Visualization Project. AT\&T Research.
  \url{http://www.graphviz.org}.
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{UG}}, \code{\link{DAG}}, \code{\link{makeMG}}, \code{\link{plotGraph}}}
\examples{
## A directed acyclic graph
d <- DAG(y1 ~ y2+y6, y2 ~ y3, y3 ~ y5+y6, y4 ~ y5+y6)
\dontrun{drawGraph(d)}

## An undirected graph
g <- UG(~giova*anto*armo + anto*arj*sara) 
\dontrun{drawGraph(d)}

## An ancestral graph
ag <- makeMG(ug=UG(~y0*y1), dg=DAG(y4~y2, y2~y1), bg=UG(~y2*y3+y3*y4))
drawGraph(ag, adjust = FALSE)
drawGraph(ag, adjust = FALSE)

## A more complex example with coordinates: the UNIX evolution
xy <-
structure(c(5, 15, 23, 25, 26, 17, 8, 6, 6, 7, 39, 33, 23, 49, 
19, 34, 13, 29, 50, 68, 70, 86, 89, 64, 81, 45, 64, 49, 64, 87, 
65, 65, 44, 37, 64, 68, 73, 85, 83, 95, 84, 0, 7, 15, 27, 44, 
37, 36, 20, 51, 65, 44, 64, 59, 73, 69, 78, 81, 90, 97, 89, 72, 
85, 74, 62, 68, 59, 52, 48, 43, 50, 34, 21, 18, 5, 1, 10, 2, 
11, 2, 1, 44), .Dim = c(41, 2), .Dimnames = list(NULL, c("x", 
"y")))
Unix <- DAG(
                SystemV.3 ~ SystemV.2,
                SystemV.2 ~ SystemV.0,
                SystemV.0 ~ TS4.0,
                TS4.0 ~ Unix.TS3.0 + Unix.TS.PP + CB.Unix.3,
                PDP11.SysV ~ CB.Unix.3,
                CB.Unix.3 ~ CB.Unix.2,
                CB.Unix.2 ~ CB.Unix.1,
                Unix.TS.PP ~ CB.Unix.3,
                Unix.TS3.0 ~ Unix.TS1.0 + PWB2.0 + USG3.0 + Interdata,
                USG3.0 ~ USG2.0,
                PWB2.0 ~ Interdata + PWB1.2,
                USG2.0 ~ USG1.0,
                CB.Unix.1 ~ USG1.0,
                PWB1.2 ~ PWB1.0,
                USG1.0 ~ PWB1.0,
                PWB1.0 ~ FifthEd,
                SixthEd ~ FifthEd,
                LSX ~ SixthEd,
                MiniUnix ~ SixthEd,
                Interdata ~ SixthEd,
                Wollongong ~ SixthEd,
                SeventhEd ~ Interdata,
                BSD1 ~ SixthEd,
                Xenix ~ SeventhEd,
                V32 ~ SeventhEd,
                Uniplus ~ SeventhEd,
                BSD3 ~ V32,
                BSD2 ~ BSD1,
                BSD4 ~ BSD3,
                BSD4.1 ~ BSD4,
                EigthEd ~ SeventhEd + BSD4.1,
                NinethEd ~ EigthEd,
                Ultrix32 ~ BSD4.2,
                BSD4.2 ~ BSD4.1,
                BSD4.3 ~ BSD4.2,
                BSD2.8 ~ BSD4.1 + BSD2,
                BSD2.9 ~ BSD2.8,
                Ultrix11 ~ BSD2.8 + V7M + SeventhEd,
                V7M ~ SeventhEd
                )
drawGraph(Unix, coor=xy, adjust=FALSE)
# dev.print(file="unix.fig", device=xfig) # Edit the graph with Xfig
}
\keyword{graphs}
\keyword{hplot}
\keyword{iplot}
