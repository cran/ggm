\name{rsphere}
\alias{rsphere}
\title{Random vectors on a sphere}
\description{
  Generates a sample of points uniformly distributed
  on the surface of a sphere in d-space.
}
\usage{
rsphere(n, d)
}
\arguments{
  \item{n}{an integer, the sample size.}
  \item{d}{an integer, the dimension of the space. For example, a circle
    is defined in 2D-space, a sphere in 3D-space.}
}
\details{The algorithm is based on normalizing to length 1 each d-vector
   of a sample from a multivariate normal \eqn{N(0, I)}.
}
\value{a matrix of \code{n} rows and \code{d} columns.
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{rnorm}}, \code{\link{rcorr}}}
\examples{
## 100 points on circle
z <- rsphere(100,2)
plot(z)

## 100 points on a sphere
z <- rsphere(100, 3)
pairs(z)
}
\keyword{distribution}
\keyword{multivariate}


