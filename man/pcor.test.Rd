\name{pcor.test}
\alias{pcor.test}
\title{Test for zero partial association}
\description{
  Test for conditional independence between two variables, given the
  other ones, assuming a multivariate normal distribution.
}
\usage{
pcor.test(r, q, n)
}
\arguments{
  \item{r}{a partial correlation coefficient, computed by \code{\link{pcor}}.}
  \item{q}{the number of variables in the conditioning set.}
  \item{n}{integer > 0, the sample size.}
}
\value{
  \item{tval}{The Student's t-test statistic.}
  \item{df}{The degrees of freedom}
  \item{pvalue}{The P-value, assuming a two-sided alternative.}
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{pcor}}, \code{\link{shipley.test}}}
\examples{
## Are 2,3 independent given 1?
data(marks)
pcor.test(pcor(c(2,3,1), var(marks)), 1, n=88)
}
\keyword{htest}
\keyword{multivariate}
