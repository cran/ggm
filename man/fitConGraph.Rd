\name{fitConGraph}
\alias{fitConGraph}
\title{Fitting a Gaussian concentration graph model}
\description{
  Fits a concentration graph (a covariance selection model).
}
\usage{
fitConGraph(amat, S, n, cli = NULL, alg = 3, pri = FALSE, tol = 1e-06)
}
\arguments{
  \item{amat}{a square Boolean matrix representing
    the adjacency matrix of an UG}
  \item{S}{the sample covariance matrix}
  \item{n}{an integer denoting the sample size}
  \item{cli}{a list containing the cliques of the graph. The components
   of the list are character vectors containing the names of the nodes in the
   cliques. The names must match the names of the vertices. The
   knowledge of the cliques is not needed. If the cliques are not specified
   the function uses the algorithm
   by Hastie et al. (2009, p. 446).}
   \item{alg}{The algorithm used.}
  \item{pri}{If TRUE is verbose}
  \item{tol}{a small positive number indicating the tolerance
    used in convergence tests.}
}
\details{
  The algorithms for fitting concentration graph
  models by maximum likelihood are discussed in Speed and Kiiveri (1986).
  If the cliques are known
  the function uses the iterative proportional
  fitting algorithm described  by Whittaker (1990, p. 184).
  If the cliques are not specified
   the function uses the algorithm
   by Hastie et al. (2009, p. 631ff).
}
\value{
\item{Shat}{the fitted covariance matrix.}
\item{dev}{the `deviance' of the model.}
\item{df}{the degrees of freedom.}
\item{it}{the iterations.}
}
\references{
    Cox, D. R. & Wermuth, N. (1996). \emph{Multivariate
      dependencies}. London: Chapman & Hall.

  Hastie, T.,  Tibshirani, R. and Friedman, J. (2009).
  \emph{The elements of statistical learning.}  Springer Verlag: New York.

  Speed, T.P. and Kiiveri, H (1986). Gaussian
  Markov distributions over finite graphs. \emph{Annals of
    Statistics}, 14, 138--150.

  Whittaker, J. (1990). \emph{Graphical models in applied
    multivariate statistics}. Chichester: Wiley.
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{UG}}, \code{\link{fitDag}},
 \code{\link{marks}}}
\examples{
## A model for the mathematics marks (Whittaker, 1990)
data(marks)
## A butterfly concentration graph
G <- UG(~ mechanics*vectors*algebra + algebra*analysis*statistics)
fitConGraph(G, cov(marks), nrow(marks))
## Using the cliques

cl = list(c("mechanics", "vectors",   "algebra"), c("algebra", "analysis" ,  "statistics"))
fitConGraph(G, S = cov(marks), n = nrow(marks), cli = cl)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
