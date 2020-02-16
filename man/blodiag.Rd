\name{blodiag}
\alias{blodiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Block diagonal matrix
}
\description{
Split a vector x into a block diagonal matrix.
}
\usage{
blodiag(x, blo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector of length \code{n}.
}
  \item{blo}{
A vector of positive integers such that \code{sum(blo) == n}.
}
}
\value{
A block-diagonal matrix with as many row as elements of \code{blo} and \code{n} columns.
The vector \code{x} is split into \code{length(blo)} sub-vectors and these are the blocks  
of the resulting matrix.
}

\author{
Giovanni M. Marchetti
}

\seealso{
\code{\link{blkdiag}}, \code{\link{diag}}
}
\examples{
blodiag(1:10, blo = c(2, 3, 5)) 
blodiag(1:10, blo = c(3,4,0,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{matrix}

