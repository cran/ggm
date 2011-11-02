\name{isADMG}
\alias{isADMG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Acyclic directed mixed graphs
}
\description{
Check if it is an adjacency matrix of an ADMG
}
\usage{
isADMG(amat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{
%%     ~~Describe \code{amat} here~~
}
}
\details{
  Checks if the following conditions must hold:
  (i) no undirected edge meets an arrowhead; (ii) no directed cycles;
}
\value{
 A logical value, \code{TRUE} if it is an ancestral graph and \code{FALSE}
otherwise.
}

\references{Richardson, T. S. and Spirtes, P. (2002). Ancestral Graph Markov
  Models. \emph{Annals of Statistics}, 30(4), 962--1030.
}
\author{
Giovanni M. Marchetti, Mathias Drton
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{makeMG}}, \code{\link{isADMG}}}   
\examples{
	## Examples from Richardson and Spirtes (2002)
	a1 <- makeMG(dg=DAG(a~b, b~d, d~c), bg=UG(~a*c))  
	isADMG(a1)    # Not an AG. (a2) p.969    
	a2 <- makeMG(dg=DAG(b ~ a, d~c), bg=UG(~a*c+c*b+b*d))           # Fig. 3 (b1) p.969  
	isADMG(a2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphs} 
\keyword{ancestral graph}      
\keyword{mixed graph}
\keyword{models}
\keyword{multivariate}
