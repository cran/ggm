\name{swp}
\alias{swp}
\title{Sweep operator}
\description{
Sweeps a covariance matrix with respect to a subset of indices. 
}
\usage{
swp(V, b)
}
\arguments{
  \item{V}{a symmetric positive definite matrix, the covariance matrix.}
  \item{b}{a subset of indices of the columns of \code{V}.}
}
\details{
The sweep operator has been introduced by Beaton (1964) as a tool for
inverting symmetric matrices (see Dempster, 1969).
}
\value{
  a square matrix \code{U} of the same order as \code{V}. If \code{a} is
  the complement of \code{b}, then \code{U[a,b]} is the matrix of
  regression coefficients of \code{a} given \code{b} and \code{U[a,a]}
  is the corresponding covariance matrix of the residuals.

  If \code{b} is empty the function returns \code{V}.

  If \code{b} is the vector \code{1:nrow(V)} (or its permutation) then
  the function returns the opposite of the inverse of \code{V}.
}
\references{
  Beaton, A.E. (1964). \emph{The use of special matrix operators
    in statistical calculus}. Ed.D. thesis, Harvard
  University. Reprinted as Educational Testing Service Research Bulletin
  64-51. Princeton.
  
  Dempster, A.P. (1969). \emph{Elements of continuous multivariate
    analysis}.  Reading: Addison-Wesley.

}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{fitDag}}}
\examples{
## A very simple example
V <- matrix(c(10, 1, 1, 2), 2, 2)
swp(V, 2)
}
\keyword{array}
\keyword{algebra}
\keyword{models}
\keyword{multivariate}
