\name{fitUg}
\alias{fitUg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gaussian Markov models specified by an UG}
\description{
  Fits a concentration graph (a covariance selection model) or
  a covariance graph to a sample covariance matrix, assuming
  a Gaussian model.
}
\usage{
fitUg(gmat, V, n, cov = FALSE, pri = FALSE, tol = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gmat}{a square Boolean matrix representing
    the edge matrix of the DAG}
  \item{V}{a symmetric positive definite matrix, the
    sample covariance matrix}
  \item{n}{an integer >0, the sample size}
  \item{cov}{a logical value. If \code{cov=TRUE}, then
    it is fitted a covariance graph. Else it is fitted
    a concentration graph.}
  \item{pri}{a logical value. If \code{TRUE} a the value
    of the criterion at each iteration is printed.}
  \item{tol}{a positive number indicating the tolerance
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
  \item{type}{a character string indicating the type of the graph,
  concentration or covariance.}
\item{cliques}{a list of numeric vectors indicating the cliques
of the UG.}
\item{Vhat}{the fitted covariance matrix.}
\item{Rhat}{the fitted correlation matrix.}
\item{Phat}{the fitted partial correlation matrix.}
\item{dev}{the `deviance' (\eqn{-2 \log L}) of the model.}
\item{df}{the degrees of freedom}
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
  \code{\link{cliques}}, \code{marks}}
\examples{
## A model for the sample covariance matrix of the
## mathematics marks (cfr. Whittaker, 1990)
data(marks)
V <- cov(marks) * 87 / 88
## A butterfly concentration graph
fitUg(UG(~ mec*vec*alg + alg*ana*sta),V , n=88)

## Correlations among four strategies to cope with stress for 
## 72 students. Cfr. Cox & Wermuth (1996), p. 73.
##  Y = cognitive avoidance
##  X = vigilance
##  V = blunting
##  U = monitoring

R <- matrix(c(
   1.00, -0.20,  0.46,  0.01,
  -0.20,  1.00,  0.00,  0.47,
   0.46,  0.00,  1.00, -0.15,
   0.01,  0.47, -0.15,  1.00), 4, 4)
nam <- c("Y", "X", "V", "U") 
dimnames(R) <- list(nam, nam)

## A chordless 4-cycle covariance graph

gr <- UG(~ Y*X + X*U + U*V + V*Y)
fitUg(gr, R, 72, cov=TRUE)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
