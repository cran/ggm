\name{anger}
\alias{anger}
\docType{data}
\title{
Anger data
}
\description{
Anger data
}
\usage{data(anger)}
\format{
A \eqn{4 \times 4} covariance matrix for the following variables.  
  \describe{
\item{Y}{}
\item{X}{}
\item{V}{}
\item{U}{}
  }   
}
\details{
See Cox and Wermuth (1996).
}
\source{
Nanny
}
\references{
  Cox, D. R. \& Wermuth, N. (1996). \emph{Multivariate
  dependencies}. London: Chapman \& Hall.}

\examples{
data(anger) 
G = UG(~ Y*X + X*V + V*U + U*Y)
fitConGraph(G,anger, 100) 
}
\keyword{datasets}
