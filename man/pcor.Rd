\name{pcor}
\alias{pcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Partial correlation}
\description{
  Computes the partial correlation between two variables given
  a set of other variables.
}
\usage{
pcor(u, S)
}
%- maybe also 'usage' for other objects documented here.
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
## The correlation between vec and alg given ana and sta
 pcor(c("vec", "alg", "ana", "sta"), var(marks))
## The same
pcor(c(2,3,4,5), var(marks))
## The correlation between vec and alg given sta
 pcor(c("vec", "alg", "sta"), var(marks))
## The marginal correlation between ana and sta 
pcor(c("ana","sta"), var(marks))
}
\keyword{models}
\keyword{multivariate}
