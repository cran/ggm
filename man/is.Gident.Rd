\name{is.Gident}
\alias{is.Gident}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{G-identifiability of an UG}
\description{
Tests if an undirected graph is G-identifiable. 
}
\usage{
is.Gident(gmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gmat}{a symmetric matrix with dimnames representing the edge
    matrix of an undirected graph}
}
\details{
  An undirected graph is G-identifiable if every connected component of
  the complementary graph contains an odd cycle (cfr. Stanghellini and
  Wermuth, 2003). The algorithm used follows suggestions by
  Tarantola and Vicard (2002).
}
\value{
  a logical value, \code{TRUE} if the graph is G-identifiable and
  \code{FALSE} if it is not.
}
\references{Stanghellini, E. \& Wermuth, N. (2003). On the
  identification of directed acyclic graph models with one
  hidden variable. Submitted and available at
  \url{http://psystat.sowi.uni-mainz.de}.
  
  Stanghellini, E. (1997). Identification of a single-factor
  model using graphical Gaussian rules. \emph{Biometrika}, 87, 199--205.

  Tarantola, C. \& Vicard, P. (2002). Spanning trees and identifiability
  of a single-factor model. \emph{Statistical Methods \& Applications},
  11, 139--152.
  
  Vicard, P. (2000). On the identification of a single-factor model
  with correlated residuals. \emph{Biometrika}, 87, 199--205.
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{UG}}, \code{\link{cycleMatrix}}}
\examples{
## A not G-identifiable UG
G1 <- UG(~ a*b + u*v)
is.Gident(G1)
## G-identifiable UG
G2 <- UG(~ a + b + u*v)
is.Gident(G2)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}

