\name{marks}
\alias{marks}
\docType{data}
\title{Mathematics marks}
\description{
Examination marks of 88 students in five subjects.
}
\usage{data(marks)}
\format{
  A data frame with 88 observations on the following 5 variables.
  \describe{
    \item{mec}{a numeric vector, mark in Mechanics}
    \item{vec}{a numeric vector, mark in Vectors}
    \item{alg}{a numeric vector, mark in Algebra}
    \item{ana}{a numeric vector, mark in Analysis}
    \item{sta}{a numeric vector, mark in Statistics }
  }
}
\details{
  Mechanics and Vectors were closed book examinations. Algebra, Analysis
  and Statistics were open book examinations.
}
\source{
  Mardia, K.V., Kent, J.T. \& Bibby, (1979). \emph{Multivariate
  analysis}. London: Academic Press.
}
\references{
  Whittaker, J. (1990). \emph{Graphical models in applied
    multivariate statistics}. Chichester: Wiley.
}
\examples{
data(marks)
pairs(marks)
}
\keyword{datasets}
