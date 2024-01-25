\name{stress}
\alias{stress}
\docType{data}
\title{
Stress
}
\description{
Stress data
}
\usage{data(stress)}
\format{
A \eqn{4 \times 4} covariance matrix for the following variables.
\describe{
\item{Y}{}
\item{V}{}
\item{X}{}
\item{U}{}
  }
}
\details{
See Cox and Wermuth (1996).
}
\references{
Cox, D. R. & Wermuth, N. (1996). \emph{Multivariate
dependencies}. London: Chapman & Hall.

Slangen K., Kleemann P.P and Krohne H.W. (1993). Coping with surgical stress. In: Krohne H. W. (ed.).  \emph{Attention and avoidance: Strategies in coping with aversiveness}. New York, Heidelberg: Springer, 321-346.
}
\examples{
data(stress)
G = UG(~ Y*X + X*V + V*U + U*Y)
fitConGraph(G, stress, 100)
}
\keyword{datasets}
