\name{fitConGraph}
\alias{fitConGraph}
\title{Fitting of Gaussian concentration graph models}
\description{
  Fits a concentration graph (a covariance selection model) 
  to a sample covariance matrix, assuming a Gaussian model.
}
\usage{
fitConGraph(amat, S, n, pri = FALSE, alg=2, tol = 1e-06)
}
\arguments{
  \item{amat}{a square Boolean matrix representing
    the adjacency matrix of the DAG}
  \item{S}{a symmetric positive definite matrix, the
    sample covariance matrix}
  \item{n}{an postive integer, the sample size}
  \item{pri}{a logical value. If \code{TRUE} a the value
    of the deviance at each iteration is printed.}
  \item{alg}{and integer specifying the algorithm: if \code{alg=1}
  then the covariance matrix is updated, if \code{alg=2} 
  the concentration matrix is updated. 
  }  
  \item{tol}{a small positive number indicating the tolerance
    used in convergence tests.
  }
}
\details{
  Algorithms for fitting Gaussian graphical models specified
  by undirected graphs are discussed in Speed \& Kiiveri (1986).
  This function is based on the iterative proportional
  fitting algorithm described on p. 184 of Whittaker (1990). 
}
\value{
\item{Shat}{the fitted covariance matrix.}
\item{dev}{the `deviance' of the model.}
\item{df}{the degrees of freedom.}
\item{it}{the iterations.}
}
\references{
    Cox, D. R. \& Wermuth, N. (1996). \emph{Multivariate
      dependencies}. London: Chapman \& Hall.
      
  Speed, T.P. \& Kiiveri, H (1986). Gaussian
  Markov distributions over finite graphs. \emph{Annals of
    Statistics}, 14, 138--150.
  
  Whittaker, J. (1990). \emph{Graphical models in applied
    multivariate statistics}. Chichester: Wiley.
}
\author{Giovanni M. Marchetti} 
\seealso{\code{\link{UG}}, \code{\link{fitDag}},
  \code{\link{cliques}}, \code{\link{marks}}}
\examples{
## A model for the sample covariance matrix of the
## mathematics marks (Whittaker, 1990)
data(marks)
S <- cov(marks) * 87 / 88
## A butterfly concentration graph
fitConGraph(UG(~ mechanics*vectors*algebra + algebra*analysis*statistics), S , n = 88)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
