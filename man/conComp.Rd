\name{conComp}
\alias{conComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Connectivity components}
\description{
Finds the connectivity components of a graph.
}
\usage{
conComp(gmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gmat}{a square Boolean matrix with dimnames, the edge matrix of
    a DAG or a UG.}
}
\value{
  an integer vector representing a partition of the set of nodes.
}
\references{Lauritzen, S. (1996). \emph{Graphical models}. Oxford:
Clarendon Press.}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{UG}}, \code{\link{anGraph}}}
\examples{
## three connected components
conComp(UG(~a*c+c*d+e+g*o*u))
## a connected graph
conComp(UG(~ a*b+b*c+c*d+d*a))
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
