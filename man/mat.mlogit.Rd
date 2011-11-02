\name{mat.mlogit}
\alias{mat.mlogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate logistic parametrization
}
\description{
 Find matrices \code{C} and \code{M} of e binary multivariate logistic 
 parameterization. 

}
\usage{
mat.mlogit(d, P = powerset(1:d))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
A positive integer, the number of binary responses.
}
  \item{P}{
A list of vectors of integers specifying margins. For instance
\code{list(1, 2, c(1,2))}. Default: the power set of \code{1:d}. 
}
}
\details{
The power set is in the order of dimensions of the sets.
}
\value{
\item{C}{A contrast matrix.}
\item{L}{A marginalization matrix.}
}
\references{
Glonek, G. J. N. and McCullagh, P. (1995). Multivariate logistic models. Journal of the Royal
Statistical Society, Ser. B 57, 533-546.
}
\author{
Giovanni M. Marchetti
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{binomial}, \link{marg.param}}
}
\examples{ 
mat.mlogit(2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{logistic model}

