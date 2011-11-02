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
Cox, D. R. \& Wermuth, N. (1996). \emph{Multivariate
dependencies}. London: Chapman \& Hall.
}
\examples{
data(stress)
G = UG(~ Y*X + X*V + V*U + U*Y)
fitConGraph(G, stress, 100)
}
\keyword{datasets}
