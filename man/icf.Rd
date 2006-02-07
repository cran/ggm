\name{icf}
\alias{icf}
\alias{icfmag}
\title{Iterative conditional fitting}
\description{
  Main algorithm for MLE fitting of Gaussian Covariance Graphs and
  Gaussian Ancestral models.
}
\usage{
icf(bi.graph, S, start = NULL, tol = 1e-06)
icfmag(mag, S, tol = 1e-06)
}
\arguments{
  \item{bi.graph}{a symmetric matrix with dimnames representing the adjacency matrix of an undirected graph.}
  \item{mag}{a square matrix representing the adjacency matrix of an
    ancestral graph (for example returned by \code{makeAG}).}
  \item{S}{a symmetric positive definite matrix, the
    sample covariance matrix. The order of the variables
  must be the same of the order of vertices in the adjacency matrix.}
  \item{start}{a symmetric matrix used as starting value
    of the algorithm. If \code{start=NULL} the starting value
    is a diagonal matrix.}
  \item{tol}{A small positive number indicating the tolerance
    used in convergence tests.}
}
\details{
  These functions are not intended to be called directly by the user.
}
\value{
  \item{Sigmahat}{the fitted covariance matrix.}
  \item{Bhat}{matrix of the fitted regression coefficients
    associated to the directed edges.}
  \item{Omegahat}{matrix of the partial covariances of the residuals
    between regression equations.}
  \item{iterations}{the number of iterations.}
}
\references{
  Drton, M. \& Richardson, T. S. (2003). A new algorithm for
  maximum likelihood estimation in Gaussian graphical models for
  marginal independence. \emph{Proceedings
    of the Ninetheen Conference on Uncertainty in Artificial
    Intelligence}, 184--191.

    Drton, M. \& Richardson, T. S. (2004). Iterative Conditional Fitting
  for Gaussian Ancestral Graph Models.  Proceedings
  of the 20th Conference on Uncertainty in Artificial Intelligence,
  Department of Statistics, 130--137.
}
\author{Mathias Drton}
\seealso{\code{\link{fitCovGraph}}, \code{\link{fitAncestralGraph}}}
\keyword{internal}
