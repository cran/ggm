\name{diagv}
\alias{diagv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix product with a diagonal matrix 
}
\description{
Computes faster the product of a diagonal matrix times a full matrix.
}
\usage{
diagv(v, M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{
A numeric vector specifying the elements on the diagonal of a matrix.
}
  \item{M}{
A numeric matrix  compatible with the product \eqn{D_v M}.
}
}
\details{
Computes \eqn{N = D_v M} where \eqn{D_v} is diagonal 
avoiding the \code{diag} operator.
}
\value{
A matrix \code{N}.
}

\seealso{
\code{\link{diag}}
}
\examples{
v <- 1:1000
M <- matrix(runif(3000), 1000, 3)
dim(diagv(v, M))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{matrix}

