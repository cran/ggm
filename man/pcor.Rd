\name{pcor}
\alias{pcor}
\title{Partial correlation}
\description{
  Computes the partial correlation between two variables given
  a set of other variables.
}
\usage{
pcor(u, S)
}
\arguments{
  \item{u}{a vector of integers of length > 1. The first two integers are
    the indices of variables the correlation of which must be
    computed. The rest of the vector is the conditioning set.
  }
  \item{S}{a symmetric positive definite matrix, a sample covariance
    matrix.}
}
\value{
  a scalar, the partial correlation matrix between variables \code{u[1]}
  and \code{u[2]} given \code{u[-c(1,2)]}.
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{cor}},  \code{\link{parcor}}, \code{\link{correlations}}}
\examples{
data(marks)
## The correlation between vectors and algebra given analysis and statistics
 pcor(c("vectors", "algebra", "analysis", "statistics"), var(marks))
## The same
pcor(c(2,3,4,5), var(marks))
## The correlation between vectors and algebra given statistics
 pcor(c("vectors", "algebra", "statistics"), var(marks))
## The marginal correlation between analysis and statistics 
pcor(c("analysis","statistics"), var(marks))
}
\keyword{models}
\keyword{multivariate}
