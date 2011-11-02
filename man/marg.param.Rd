\name{marg.param}
\alias{marg.param}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Link function of marginal log-linear parameterization
}
\description{
Provides the contrast and marginalization matrices for the marginal parametrization
of a probability vector.
}
\usage{
marg.param(lev, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lev}{
Integer vector containing the number of levels of each variable.
}
  \item{type}{
A character vector   with elements \code{"l"}, \code{"g"}, \code{"c"}, or 
\code{"r"} indicating the type of logit. The meaning is as follows:
\code{"g"} for global, \code{"c"} for continuation,
\code{"r"} for reverse continuation and \code{"l"} for local.
}
}
\details{
See Bartolucci, Colombi and Forcina (2007).
}
\value{
\item{C}{Matrix of constrasts (the first \code{sum(lev)-length(r)} elements are
      referred to univariate logits)}
  \item{M}{Marginalization matrix with elements 0 and 1.}
\item{G}{Corresponding design matrix for the corresponding log-linear model.}
}
\references{
Bartolucci, F., Colombi, R. and Forcina, A. (2007). An extended class 
of marginal link functions for modelling contingency tables 
by equality and inequality constraints. Statist. Sinica 17, 691-711.
}
\author{
Francesco Bartolucci, Antonio Forcina, Giovanni M. Marchetti
}
\note{
Assumes that the vector of probabilities is in inv lex order. 
The interactions are returned in order of dimension, like e.g., 
 1, 2, 3, 12, 13, 23, 123.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mat.mlogit}}
}
\examples{    
marg.param(c(3,3), c("l", "g"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{logistic models}
\keyword{ordinal models}% __ONLY ONE__ keyword per line
