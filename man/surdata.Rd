\name{surdata}
\alias{surdata}
\docType{data}
\title{
A simulated data set 
}
\description{
Simulated data following a seemingly unrelated regression model.  
}
\usage{data(surdata)}
\format{
  A data frame with 600 observations on the following 4 variables.
  \describe{
    \item{\code{A}}{a numeric response vector}
    \item{\code{B}}{a numeric response vector}
    \item{\code{X}}{a numeric vector}
    \item{\code{Z}}{a numeric vector with codes \code{1} and \code{-1} for a binary variables.}
  }
}

\examples{
data(surdata)
pairs(surdata)
}
\keyword{datasets}
