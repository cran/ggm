\name{unmakeAG}
\alias{unmakeAG}
\title{Ancestral Graph components}
\description{
  Splits the adjacency matrix of an ancestral graph into
  undirected and bidirected components.
}
\usage{
unmakeAG(mag)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mag}{a square matrix, with dimnames, representing an
  ancestral graph.}
}
\value{
  It is the inverse of \code{makeAG}. It returns the following
  components.
  \item{dag}{the adjacency matrix of the directed edges.}
  \item{ug}{the adjacency matrix of the undirected edges.}
  \item{bg}{the adjacency matrix of the bidirected edges.}
}
\author{Mathias Drton}
\seealso{\code{\link{makeAG}}}
\examples{
ag <- makeAG(ug=UG(~y0*y1), dag=DAG(y4~y2, y2~y1), bg=UG(~y2*y3+y3*y4))
unmakeAG(ag)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
