\name{isGident}
\alias{isGident}
\title{G-identifiability of an UG}
\description{
Tests if an undirected graph is G-identifiable. 
}
\usage{
isGident(amat)
}
\arguments{
  \item{amat}{a symmetric matrix with dimnames representing the adjacency
    matrix of an undirected graph}
}
\details{
  An undirected graph is said G-identifiable if every connected component of
  the complementary graph contains an odd cycle (Stanghellini and
  Wermuth, 2005). See also Tarantola and Vicard (2002).
}
\value{
  a logical value, \code{TRUE} if the graph is G-identifiable and
  \code{FALSE} if it is not.
}
\references{
  Stanghellini, E. \& Wermuth, N. (2005). On the
  identification of path-analysis models with one
  hidden variable. \emph{Biometrika}, 92(2), 337-350.

  Stanghellini, E. (1997). Identification of a single-factor
  model using graphical Gaussian rules. \emph{Biometrika}, 84, 241--244.

  Tarantola, C. \& Vicard, P. (2002). Spanning trees and identifiability
  of a single-factor model. \emph{Statistical Methods \& Applications},
  11, 139--152.
  
  Vicard, P. (2000). On the identification of a single-factor model
  with correlated residuals. \emph{Biometrika}, 87, 199--205.
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{UG}}, \code{\link{cmpGraph}}, \code{\link{cycleMatrix}}}
\examples{
## A not G-identifiable UG
G1 <- UG(~ a*b + u*v)
isGident(G1)
## G-identifiable UG
G2 <- UG(~ a + b + u*v)
isGident(G2)
## G-identifiable UG
G3 <- cmpGraph(UG(~a*b*c+x*y*z))
isGident(G3)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}

