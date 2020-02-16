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
    \item{mechanics}{a numeric vector, mark in Mechanics}
    \item{vectors}{a numeric vector, mark in Vectors}
    \item{algebra}{a numeric vector, mark in Algebra}
    \item{analysis}{a numeric vector, mark in Analysis}
    \item{statistics}{a numeric vector, mark in Statistics }
  }
}
\details{
  Mechanics and Vectors were closed book examinations. Algebra, Analysis
  and Statistics were open book examinations.
}
\source{
  Mardia, K.V., Kent, J.T. and Bibby, (1979). \emph{Multivariate
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
