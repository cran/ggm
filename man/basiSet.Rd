\name{basiSet}
\alias{basiSet}
\title{Basis set of a DAG}
\description{
  Finds a basis set for the conditional independencies
  implied by a directed acyclic graph, that is a minimal set of
  independencies that imply all the other ones.
}
\usage{
basiSet(amat)
}
\arguments{
  \item{amat}{a square matrix with dimnames representing the adjacency matrix
    of a DAG.}
}
\details{
  Given a DAG and a pair of non adjacent nodes
  \eqn{(i,j)} such that \eqn{j} has higher causal order than \eqn{i},
  the set of independency statements \eqn{i} independent of
  \eqn{j} given the union of the parents of both \eqn{i} and \eqn{j}
  is a basis set (see Shipley, 2000). This basis set has the property
  to lead to independent test statistics.
}
\value{
  a list of vectors representing several conditional independence
  statements. Each vector contains the names of two non adjacent
  nodes followed by the names of nodes in the conditioning set (which
  may be empty).
}
\references{Shipley, B. (2000). A new inferential test
for path models based on directed acyclic graphs. \emph{Structural
  Equation Modeling}, 7(2), 206--218.}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{shipley.test}}, \code{\link{dSep}}, \code{\link{DAG}}}
\examples{
## See Shipley (2000), Figure 2, p. 213
A <- DAG(x5~ x3+x4, x3~ x2, x4~x2, x2~ x1)
basiSet(A)
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}
