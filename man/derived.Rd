\name{derived}
\alias{derived}
\docType{data}
\title{Data on blood pressure body mass and age
}
\description{
Raw data on blood pressure, body mass and age on 44 
female patients, and covariance matrix for derived variables.
}
\usage{data(derived)}
\format{A list containing a dataframe
\code{raw} with 44 lines and 5 columns and a
 symmetric 4x4 covariance matrix \code{S}.
 
 The following is the description of the variables
 in the dataframe \code{raw}
  \describe{
    \item{\code{Sys}}{Systolic blood pressure, in mm Hg}
    \item{\code{Dia}}{Diastolic blood pressure, in mm Hg}
    \item{\code{Age}}{Age of the patient, in years}
    \item{\code{Hei}}{Height, in cm}
    \item{\code{Wei}}{Weight, in kg}  
  }
The following is the description of the variables
for the covariance matrix \code{S}.  
\describe{
\item{\code{Y}}{Derived variable \code{Y=log(Sys/Dia)}}
\item{\code{X}}{Derived variables \code{X=log(Dia)}}
\item{\code{Z}}{Body mass index \code{Z=Wei/(Hei/100)^2}}
\item{\code{W}}{Age}
}
}
\references{
  Wermuth N. and Cox D.R. (1995). Derived 
  variables calculated from similar joint responses:
  some characteristics and examples. \emph{Computational Statistics and Data
  Analysis}, 19, 223-234.
}
\examples{
# A DAG model with a latent variable U
G = DAG(Y ~ Z + U, X ~ U + W, Z ~ W)

data(derived)

# The model fitted using the derived variables
out = fitDagLatent(G, derived$S, n = 44, latent = "U")

# An ancestral graph model marginalizing over U
H = AG(G, M = "U")

# The ancestral graph model fitted obtaining the 
# same result
out2 = fitAncestralGraph(H, derived$S, n = 44)
}
\keyword{datasets}
