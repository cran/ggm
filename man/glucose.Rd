\name{glucose}
\alias{glucose}
\docType{data}
\title{Glucose control}
\description{
Data on glucose control of diabetes patients.
}
\usage{data(glucose)}
\format{
  A data frame with 68 observations on the following 8 variables.
  \describe{
    \item{Y}{a numeric vector, Glucose control (glycosylated
    haemoglobin), values up to about 7 or 8 indicate good glucose control.}
    \item{X}{a numeric vector, a score for knowledge about the illness.}
    \item{Z}{a numeric vector, a score for fatalistic externality (mere
      chance determines what occurs).}
    \item{U}{a numeric vector, a score for social externality (powerful
      others are responsible).}
    \item{V}{a numeric vector, a score for internality (the patient is
    him or herself responsible).}
    \item{W}{a numeric vector, duration of the illness in years.}
    \item{A}{a numeric vector, level of education, with levels \code{-1}: at least
      13 years of formal schooling, \code{1}: less then 13 years.}
    \item{B}{a numeric vector, gender with levels \code{-1}: females, \code{1}: males.}
  }
}
\details{
  Data on 68 patients with fewer than 25 years of diabetes. They were
  collected at the University of Mainz to identify psychological
  and socio-economic variables possibly important for glucose control,
  when patients choose the appropriate dose of treatment depending
  on the level of blood glucose measured several times per day.

  The variable   of primary interest is \code{Y}, glucose control, measured
  by  glycosylated haemoglobin. \code{X}, knowledge about the illness,
  is a response of secondary interest. Variables \code{Z}, \code{U} and
  \code{V} measure patients' type of attribution, called fatalistic
  externality, social externality and internality. These are intermediate
  variables. Background   variables are \code{W}, the duration of the
  illness, \code{A} the duration of formal schooling and \code{B},
  gender. The background variables \code{A} and \code{B} are binary
  variables with coding \code{-1}, \code{1}. 
}
\source{
Cox & Wermuth (1996), p. 229.
}
\references{Cox, D. R. \& Wermuth, N. (1996). \emph{Multivariate
  dependencies}. London: Chapman \& Hall.}
\examples{
data(glucose)
## See Cox & Wermuth (1996), Figure 6.3 p. 140
coplot(Y ~ W | A, data=glucose)
}
\keyword{datasets}
