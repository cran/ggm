\name{conComp}
\alias{conComp}
\title{Connectivity components}
\description{
Finds the connectivity components of a graph.
}
\usage{
conComp(amat, method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{a square matrix with dimnames, the adjacency matrix of an UG.}
  \item{method}{an integer 1 or 2 to choose the method used to find the components. Method 2 is more efficient for large graphs.}
}
\value{
  an integer vector representing a partition of the set of nodes.
}
\references{Lauritzen, S. (1996). \emph{Graphical models}. Oxford:
Clarendon Press.}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{UG}}}
\examples{
## three connected components
conComp(UG(~a*c+c*d+e+g*o*u))
## a connected graph
conComp(UG(~ a*b+b*c+c*d+d*a))
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
