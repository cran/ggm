\name{powerset}
\alias{powerset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power set
}
\description{
Finds the list of all subsets of a set.
}
\usage{
powerset(set, sort = TRUE, nonempty = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{set}{
A numeric or character vector.
}
  \item{sort}{
Logical value. If \code{TRUE} the subsets are sorted according 
to dimension.   Default is \code{TRUE}. 
}
  \item{nonempty}{
 Logical value. If \code{TRUE} the empty set is omitted.   Default is \code{TRUE}.
}
}
\details{
If \code{sort == FALSE} the sets are in inverse lexicographical order.  
}
\value{
A list of all subsets of \code{set}.
}

\author{
Giovanni M. Marchetti
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
powerset(c("A", "B", "C"), nonempty = FALSE)  
powerset(1:3, sort = FALSE, nonempty = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sets}

