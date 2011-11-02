\name{blkdiag}
\alias{blkdiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Block diagonal matrix
}
\description{
Block diagonal concatenation of input arguments.
}
\usage{
blkdiag(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
Variable number of matrices \code{M1, M2, ...}.
}
}
\value{
A block diagonal matrix \code{diag(M1, M2, ...)}.
}

\author{
Giovanni M. Marchetti
}

\seealso{
\code{\link{diag}}
}
\examples{
X <- c(1,1,2,2); Z <- c(10, 20, 30, 40); A <- factor(c(1,2,2,2))
blkdiag(model.matrix(~X+Z), model.matrix(~A))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{matrix}

