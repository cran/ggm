\name{fitAncestralGraph}
\alias{fitAncestralGraph}
\title{Fitting of Gaussian Ancestral Graph Models}
\description{
Iterative conditional fitting of Gaussian Ancestral Graph Models.
}
\usage{
fitAncestralGraph(amat, S, n, tol = 1e-06)
}
\arguments{
  \item{amat}{a square matrix, representing the adjacency matrix of
    an ancestral graph. }
  \item{S}{a symmetric positive definite matrix with dimnames, the
    sample covariance matrix.}
  \item{n}{the sample size, a positive integer.}
  \item{tol}{a small positive number indicating the tolerance
    used in convergence checks.}
}
\details{
  Ancestral graph models were introduced by Richardson \&
  Spirtes (2002) as a class of graphical models
  whose global Markov property is closed under conditioning and
  marginalization. In the Gaussian case, the models can be parameterized
  using precision parameters, regression coefficients, and error
  covariances (compare Richardson \& Spirtes, 2002, Section 8). This
  function finds the MLE \eqn{\hat \Lambda}{\code{Lhat}} of the precision 
  parameters by fitting a concentration  
  graph model. The MLE \eqn{\hat B}{\code{Bhat}} of the regression coefficients and
  the MLE \eqn{\hat\Omega}{\code{Ohat}} of the error covariances are obtained by
  iterative conditional fitting (Drton \& Richardson, 2003a, b). The
  three sets of parameters are 
  combined to the MLE \eqn{\hat\Sigma}{\code{Shat}} of the covariance matrix by
  matrix multiplication: 
  \deqn{\hat\Sigma = \hat B^{-1}(\hat \Lambda+\hat\Omega)\hat
    B^{-T}.}{\code{Shat = solve(Bhat) \%*\% (Lhat+Ohat) \%*\% t(solve(Bhat))}.}
  Note that in Richardson \& Spirtes (2002), the matrices \eqn{\Lambda}{Lhat}
  and \eqn{\Omega}{Ohat} are defined as submatrices.
}
\value{
  \item{Shat}{the fitted covariance matrix.}
  \item{Lhat}{matrix of the fitted precisions associated with undirected
    edges and vertices that do not have an arrowhead pointing at them.}
  \item{Bhat}{matrix of the fitted regression coefficients
    associated to the directed edges.  Precisely said \code{Bhat}
    contains ones on the diagonal and the off-diagonal entry
    \eqn{(i,j)}{(i,j)} equals the \emph{negated} MLE of the regression
    coefficient for variable \eqn{j}{j} in the regression of variable
    \eqn{i}{i} on its parents. Note that this \eqn{(i,j)}{(i,j)} entry
    in \code{Bhat} corresponds to a directed edge \eqn{j \to i}{j -> i},
    and thus to a one as \eqn{(j,i)}{(j,i)} entry in the adjacency matrix.}
  \item{Ohat}{matrix of the error covariances and variances of the residuals 
    between regression equations associated with bidirected edges and
    vertices with an arrowhead pointing at them.}
  \item{dev}{the `deviance' of the model.}
  \item{df}{the degrees of freedom.}
  \item{it}{the iterations.}
}
\references{
  Drton, M. \& Richardson, T. S. (2003a). Iterative Conditional Fitting
  for Gaussian Ancestral Graph Models.  Department of Statistics,
  University of Washington, Technical Report 437, under preparation.
  (Compare also \url{http://www.math.auc.dk/gr/gr2003/material/drton.pdf})
 
  Drton, M. \& Richardson, T. S. (2003b). A new algorithm for
  maximum likelihood estimation in Gaussian graphical models for
  marginal independence. \emph{Proceedings
    of the Nineteenth Conference on Uncertainty in Artificial
    Intelligence}, 184--191.

  Richardson, T. Spirtes, P. (2002). Ancestral Graph Markov
  Models. \emph{Annals of Statistics}. 30, 4, 962--1030.
}
\author{Mathias Drton}
\seealso{\code{\link{fitCovGraph}}, \code{\link{icf}},
  \code{\link{makeAG}}, \code{\link{fitDag}}
  }
\examples{
## A covariance matrix
"S" <- structure(c(2.93, -1.7, 0.76, -0.06,
                  -1.7, 1.64, -0.78, 0.1,
                   0.76, -0.78, 1.66, -0.78,
                  -0.06, 0.1, -0.78, 0.81), .Dim = c(4,4),
                 .Dimnames = list(c("y", "x", "z", "u"), c("y", "x", "z", "u")))
## The following should give the same fit.   
## Fit an ancestral graph y -> x <-> z <- u
fitAncestralGraph(ag1 <- makeAG(dag=DAG(x~y,z~u), bg = UG(~x*z)), S, n=100)

## Fit an ancestral graph y <-> x <-> z <-> u
fitAncestralGraph(ag2 <- makeAG(bg= UG(~y*x+x*z+z*u)), S, n=100)

## Fit the same graph with fitCovGraph
fitCovGraph(ag2, S, n=100)    

## Another example for the mathematics marks data

data(marks)
S <- var(marks)
mag1 <- makeAG(bg=UG(~mechanics*vectors*algebra+algebra*analysis*statistics))
fitAncestralGraph(mag1, S, n=88)

mag2 <- makeAG(ug=UG(~mechanics*vectors+analysis*statistics),
               dag=DAG(algebra~mechanics+vectors+analysis+statistics))
fitAncestralGraph(mag2, S, n=88) # Same fit as above
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}







