\name{checkIdent}
\alias{checkIdent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identifiability of a model with one latent variable}
\description{
  Checks  four sufficient conditions for identifiability of a Gaussian DAG model
with one latent variable.
}
\usage{
checkIdent(gmat, hidden)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gmat}{a square Boolean matrix with dimnames, the edge matrix of
    a DAG.}
  \item{hidden}{an integer representing the latent variables among the
    nodes, or the name of the node.}
}
\details{
  Stanghellini and Wermuth (2003) give some
  sufficient conditions for checking if a Gaussian
  model that factorizes according to a DAG is identified
  when there is one hidden node over which we marginalize.
  Specifically, the function checks the conditions of
  Theorem 1, (i) and (ii) and of Theorem 2 (i) and (ii).
}
\value{
  a vector of logical values, indicating if each of the  four 
  conditions of theorems 1 and 2 in Stanghellini \& Wermuth are
  \code{TRUE} or \code{FALSE}.
}
\references{
  Stanghellini, E. \& Wermuth, N. (2003). On the
  identification of directed acyclic graph models with one
  hidden variable. Submitted and available at
  \url{http://psystat.sowi.uni-mainz.de}.
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{is.Gident}},  \code{\link{InducedGraphs}}}
\examples{
## See DAG in Figure 4 (a) in Stanghellini & Wermuth (2003)
d <- DAG(y1 ~ y3, y2 ~ y3 + y5, y3 ~ y4 + y5, y4 ~ y6)
checkIdent(d, "y3")  # Identifiable
checkIdent(d, "y4")  # Not identifiable

## See DAG in Figure 5 (a) in Stanghellini & Wermuth (2003)
d <- DAG(y1 ~ y5+y4, y2 ~ y5+y4, y3 ~ y5+y4)
checkIdent(d, "y4")  # Identifiable
checkIdent(d, "y5")  # Identifiable

## A simple function to check identifiability for each node
"is.ident" <-
function(gmat, compact=TRUE){
### Check suff. conditions marginalising on each node of a DAG. 
    a <- rownames(gmat)
    for(i in a) {
        b <- checkIdent(gmat, hidden=i)
      if(compact)
        b <- any(b)
      
      cat("Node: ", i, "\n")
      print(b)
    }
  }

}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}

