\name{fitDag}
\alias{fitDag}
\title{Fitting of Gaussian DAG models}
\description{
  Fits linear recursive regressions with independent residuals specified
  by a DAG.
}
\usage{
fitDag(amat, S, n)
}
\arguments{
  \item{amat}{a square matrix with dimnames representing
    the adjacency matrix of the DAG}
  \item{S}{a symmetric positive definite matrix, the
    sample covariance matrix}
  \item{n}{an integer > 0, the sample size}
}
\details{
  \code{fitDag} checks if the order of the nodes in adjacency matrix
  is the same of \code{S} and if not it reorders the adjacency matrix
  to match the order of the variables in \code{S}. The nodes
  of the adjacency matrix may form a subset of the variables in \code{S}.
}
\value{
  \item{Shat}{the fitted covariance matrix.}
  \item{Ahat}{a square matrix of the fitted regression coefficients. The entry
    \code{Ahat[i,j]} is minus the regression coefficient of variable
    \code{i} in the regression equation \code{j}. Thus there is
    a non zero partial regression coefficient \code{Ahat[i,j]}
    corresponding to each non zero value  \code{amat[j,i]} in
    the adjacency matrix.}
  \item{Dhat}{a vector containing the partial variances of each
    variable given the parents.}
  \item{dev}{the `deviance' of the model.}
  \item{df}{the degrees of freedom.}
}
\references{Cox, D. R. & Wermuth, N. (1996). \emph{Multivariate
  dependencies}. London: Chapman & Hall.}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{DAG}}, \code{\link{swp}}.}
\examples{
dag <- DAG(y ~ x+u, x ~ z, z ~ u)
"S" <- structure(c(2.93, -1.7, 0.76, -0.06,
                   -1.7, 1.64, -0.78, 0.1,
                    0.76, -0.78, 1.66, -0.78,
                    -0.06, 0.1, -0.78, 0.81), .Dim = c(4,4),
         .Dimnames = list(c("y", "x", "z", "u"), c("y", "x", "z", "u")))
fitDag(dag, S, 200)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
