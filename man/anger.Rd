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
A covariance matrix for 4 variables measured on 684 female students. 
  \describe{
\item{X}{anxiety state}
\item{Y}{anger state}
\item{Z}{anxiety trait}
\item{U}{anger trait}
  }   
}
\details{ Trait variables are viewed as stable personality characteristics, and
state variables denote behaviour in specific situations. 
See Cox and Wermuth (1996).
}
\references{
  Cox, D. R. and Wermuth, N. (1996). \emph{Multivariate
  dependencies}. London: Chapman and Hall.
  
  Cox, D.R. and Wermuth, N. (1990). \emph{An approximation to maximum likelihood
estimates in reduced models}. 77(4), 747-761. 
}

\examples{ 
# Fit a chordless 4-cycle model 
data(anger) 
G = UG(~ Y*X + X*Z + Z*U + U*Y)
fitConGraph(G,anger, 684) 
}
\keyword{datasets}
