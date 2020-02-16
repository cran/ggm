\name{fitmlogit}
\alias{fitmlogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate logistic models
}
\description{
Fits a logistic regression model to multivariate binary responses.
}
\usage{
fitmlogit(..., C = c(), D = c(), data, mit = 100, ep = 1e-80, acc = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
Model formulae of marginal logistic models for each response and for each 
association parameters (log-odds ratios).
}
  \item{C}{
Matrix of equality constraints. 
}
  \item{D}{
Matrix of inequality cosntraints.
}
  \item{data}{
A data frame containing the responses and the explanatory variables.
}
  \item{mit}{
A positive integer: maximum number of iterations.  Default: \code{100}.
}
  \item{ep}{
A tolerance used in the algorithm: default \code{1e-80}. 
}
  \item{acc}{
A tolerance used in the algorithm: default \code{1e-4}.      
}
}
\details{
See Evans and Forcina (2011). 
}
\value{
\item{LL}{The maximized log-likelihood.}
\item{be}{The vector of the Maximum likelihood estimates of the parameters.}  
\item{S}{The estimated asymptotic covariance matrix.}
\item{P}{The estimated cell probabilities for each individual.} 
}
\references{
Evans, R.J. and Forcina, A. (2013). Two algorithms for fitting constrained marginal models. \emph{Computational Statistics and Data Analysis}, 66, 1-7. 
}
\author{
Antonio Forcina, Giovanni M. Marchetti
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{glm}}
}
\examples{    
data(surdata)                     
out1 <- fitmlogit(A ~X, B ~ Z, cbind(A, B) ~ X*Z, data = surdata)     
out1$beta
out2 <- fitmlogit(A ~X, B ~ Z, cbind(A, B) ~ 1, data = surdata)        
out2$beta
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{logistic model}% __ONLY ONE__ keyword per line
