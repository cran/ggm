\name{fitDag}
\alias{fitDag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gaussian Markov models specified by a DAG}
\description{
  Fits linear recursive regressions with independent residuals specified
  by a DAG.
}
\usage{
fitDag(gmat, Syy, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gmat}{a square Boolean matrix representing
    the edge matrix of the DAG}
  \item{Syy}{a symmetric positive definite matrix, the
    sample covariance matrix}
  \item{n}{an integer >0, the sample size}
}
\details{
  \code{fitDag} checks if the order of the nodes in edge matrix
  is the same of \code{Syy} and if not it reorders the edge matrix
  to match the order of the variables in \code{Syy}. The nodes
  of the edge matrix may form a subset of the variables in \code{Syy}.
}
\value{
  \item{A}{the matrix of the fitted regression coefficients. Each
    coefficient has the sign reversed.}
  \item{B}{the inverse of matrix B.}
  \item{Delta}{a vector containing the partial covariances of each
    variable given the parents.}
  \item{Shat}{the fitted covariance matrix.}
  \item{Khat}{the fitted concentration matrix.}
  \item{dev}{the `deviance' (\eqn{-2 \log L}) of the model.}
  \item{df}{the degrees of freedom.}
}
\references{Cox, D. R. \& Wermuth, N. (1996). \emph{Multivariate
  dependencies}. London: Chapman \& Hall.}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{DAG}}, \code{\link{swp}}.}
\examples{
dag <- DAG(y ~ x+u, x ~ z, z ~ u)
"V" <- structure(c(2.93, -1.7, 0.76, -0.06,
                   -1.7, 1.64, -0.78, 0.1,
                    0.76, -0.78, 1.66, -0.78,
                    -0.06, 0.1, -0.78, 0.81), .Dim = c(4,4),
         .Dimnames = list(c("y", "x", "z", "u"), c("y", "x", "z", "u")))
fitDag(dag, V, 200)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
