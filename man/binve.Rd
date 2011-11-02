\name{binve}
\alias{binve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inverts a marginal log-linear parametrization
}
\description{
Inverts a marginal log-linear parametrization.   
}
\usage{
binve(eta, C, M, G, maxit = 500, print = FALSE, tol = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta}{
a vector of dimension \code{t-1} where \code{t} is the number of cells
of a contingency table.  
}
  \item{C}{
A contrast matrix. 
}
  \item{M}{
A marginalization matrix.
}
  \item{G}{
G is the model matrix of the loglinear parameterization with no constant term.
}
  \item{maxit}{
an integer, specifying the maximum number of iterations. Default 500. 
}
  \item{print}{
a logical value: if \code{TRUE}, prints the criterion after each cycle.  
}
  \item{tol}{
A small value specifying the tolerance for the convergence criterion. Default:
\code{1e-10}.
}
}
\details{
A marginal log-linear link is defined by \eqn{\eta = C  (M \log p)}. See 
Bartolucci et al. (2007).
}
\value{
A vector of probabilities \code{p}.
}
\references{
Bartolucci, F., Colombi, R. and Forcina, A. (2007). An extended class 
of marginal link functions for modelling contingency tables 
by equality and inequality constraints. Statist. Sinica 17, 691-711.
}
\author{
Antonio Forcina, Giovanni M. Marchetti
}
\note{
From a Matlab function by A. Forcina, University of Perugia, Italy.
}

\seealso{
\code{\link{mat.mlogit}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{marginal log-linear models}
\keyword{discrete data}
