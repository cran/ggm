\name{fitDagLatent}
\alias{fitDagLatent}
\title{Fitting Gaussian DAG models with one latent variable}
\description{
Fits by maximum likelihood a Gaussian DAG model where one of the nodes
of the graph is latent and it is marginalised over.
}
\usage{
fitDagLatent(amat, Syy, n, latent, norm = 1, seed,
             maxit = 9000, tol = 1e-06, pri = FALSE)
}
\arguments{
  \item{amat}{a square matrix with dimnames representing
    the adjacency matrix of the DAG.}
  \item{Syy}{a symmetric positive definite matrix,
    with dimnames, the sample covariance
    matrix of the observed variables.
    The set of the observed nodes of the graph must be a subset
    of the set of the names of the variables in \code{Syy}.}
  \item{n}{a positive integer, the sample size.}
  \item{latent}{the name of the latent variable.}
  \item{norm}{an integer, the kind of normalization of the latent
    variable.
    If \code{norm=1}, the latent is scaled to have unit variance. If
    \code{norm=2}, the latent is scaled to have unit partial variance
    given its parents.}
  \item{seed}{an integer, used by \code{set.seed} to specify a random
  starting point of the EM algorithm.}
\item{maxit}{an integer denoting the maximum number of iterations
  allowed for the EM algorithm. If the convergence criterion is not
  satisfied within maxit iterations the algorithms stops and a warning
  message is returned.}
\item{tol}{a small real value, denoting the tolerance used in testing convergence.}
\item{pri}{logical, if \code{pri=TRUE} then the value of the deviance at
  each iteration is printed.}
}
\details{
For the EM algorithm used see Kiiveri (1987).  
}
\value{
  \item{Shat}{
    the fitted covariance matrix of all the variables
    including the latent one. The latent variable is the last.
    If \code{norm=1} then the variance of the latent variable
    is constrained to 1.
  }
  \item{Ahat}{
   a square matrix of the fitted regression coefficients. The entry
    \code{Ahat[i,j]} is minus the regression coefficient of variable
    \code{i} in the regression equation \code{j}. Thus there is
    a non zero partial regression coefficient \code{Ahat[i,j]}
    corresponding to each non zero value  \code{amat[j,i]} in
    the adjacency matrix.
  }
  \item{Dhat}{
    a vector containing the partial variances of each
    variable given the parents. If \code{norm=2} then the
    partial variance of the latent variable is constrained to 1. 
  }
  \item{dev}{
    the `deviance' of the model.
  }
  \item{df}{
    the degrees of freedom of the model.
  }
  \item{it}{
    the number of EM algorithm iterations at convergence.
  }
}
\references{
  Kiiveri,H. T. (1987). An incomplete data approach to the
  analysis of covariance structures. \emph{Psychometrika}, 52, 4,
  539--554.

  Joreskog, K.G. and Goldberger, A.S. (1975). Estimation of a model
  with multiple indicators and multiple causes of a single latent
  variable.
  \emph{Journal of the American Statistical Association}, 10, 631--639.
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{fitDag}}, \code{\link{checkIdent}}}
\examples{
## data from Joreskog and Goldberger (1975)
V <- matrix(c(1,     0.36,   0.21,  0.10,  0.156, 0.158,
              0.36,  1,      0.265, 0.284, 0.192, 0.324,
              0.210, 0.265,  1,     0.176, 0.136, 0.226,
              0.1,   0.284,  0.176, 1,     0.304, 0.305, 
              0.156, 0.192,  0.136, 0.304, 1,     0.344,
              0.158, 0.324,  0.226, 0.305, 0.344, 1),     6,6)
nod <- c("y1", "y2", "y3", "x1", "x2", "x3")
dimnames(V) <- list(nod,nod)
dag <- DAG(y1 ~ z, y2 ~ z, y3 ~ z, z ~ x1 + x2 + x3, x1~x2+x3, x2~x3) 
fitDagLatent(dag, V, n=530, latent="z", seed=4564)
fitDagLatent(dag, V, n=530, latent="z", norm=2, seed=145)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
