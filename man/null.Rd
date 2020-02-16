\name{null}
\alias{null}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Null space of a matrix
}
\description{
Given a matrix \code{M} find a matrix \code{N}  such that  
\eqn{N^T M} is zero. 
}
\usage{
null(M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
A matrix.
}
}

\value{
The matrix \code{N} with the basis for the null space, or an empty vector if 
the matrix \code{M} is square and of maximal rank.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Null}}, ~~~
}
\examples{
 null(c(1,1,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{matrix}

