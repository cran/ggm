\name{rcorr}
\alias{rcorr}
\title{Random correlation matrix}
\description{
Generates a random correlation matrix 
with the method of Marsaglia and Olkin (1984). 
}
\usage{
rcorr(d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{an integer > 0,  the order of the correlation matrix.}
}
\details{
  The algorithm uses \code{\link{rsphere}} to generate \eqn{d}
  vectors on a sphere in  \eqn{d}-space. If \eqn{Z}
  is a matrix with such vectors as rows, then the random
  correlation matrix is \eqn{ZZ'}.
}
\value{
a correlation matrix of order \code{d}.
}
\references{Marshall, G.\& Olkin, I. (1984).Generating correlation
  matrices. \emph{SIAM J. Sci. Stat. Comput.}, 5, 2, 470--475.}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{rsphere}}}
\examples{
## A random correlation matrix of order 3
rcorr(3)
## A random correlation matrix of order 5
rcorr(5)
}
\keyword{distribution}
\keyword{multivariate}
