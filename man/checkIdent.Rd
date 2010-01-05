\name{checkIdent}
\alias{checkIdent}
\title{Identifiability of a model with one latent variable}
\description{
  Checks  four sufficient conditions for
  identifiability of a Gaussian DAG model
  with one latent variable.
}
\usage{
checkIdent(amat, latent)
}
\arguments{
  \item{amat}{a square matrix with dimnames,
    representing the adjacency matrix of a DAG.}
  \item{latent}{an integer representing the latent variables among the
    nodes, or the name of the node.}
}
\details{
  Stanghellini and Wermuth (2005) give some
  sufficient conditions for checking if a Gaussian
  model that factorizes according to a DAG is identified
  when there is one hidden node over which we marginalize.
  Specifically, the function checks the conditions of
  Theorem 1, (i) and (ii) and of Theorem 2 (i) and (ii).
}
\value{
  a vector of length  four, indicating if the model is identified
  according to the conditions of theorems 1 and 2 in Stanghellini
  \& Wermuth (2005). The answer is \code{TRUE} if the  condition holds
  and thus the model is globally identified or \code{FALSE} if the
  condition fails, and thus we do not know if the model is identifiable.
}
\references{
  Stanghellini, E. \& Wermuth, N. (2005). On the
  identification of path-analysis models with one
  hidden variable. \emph{Biometrika}, 92(2), 337-350.
}
\author{Giovanni M. Marchetti}
\seealso{\code{\link{isGident}},  \code{\link{InducedGraphs}}}
\examples{
## See DAG in Figure 4 (a) in Stanghellini & Wermuth (2005)
d <- DAG(y1 ~ y3, y2 ~ y3 + y5, y3 ~ y4 + y5, y4 ~ y6)
checkIdent(d, "y3")  # Identifiable
checkIdent(d, "y4")  # Not identifiable?

## See DAG in Figure 5 (a) in Stanghellini & Wermuth (2005)
d <- DAG(y1 ~ y5+y4, y2 ~ y5+y4, y3 ~ y5+y4)
checkIdent(d, "y4")  # Identifiable
checkIdent(d, "y5")  # Identifiable

## A simple function to check identifiability for each node

is.ident <- function(amat){
### Check suff. conditions on each node of a DAG.
   p <- nrow(amat)
   ## Degrees of freedom
     df <- p*(p+1)/2 - p  - sum(amat==1) - p + 1
   if(df <= 0)
       warning(paste("The degrees of freedom are ", df))
    a <- rownames(amat)
    for(i in a) {
      b <- checkIdent(amat, latent=i)
      if(TRUE \%in\% b)
        cat("Node", i, names(b)[!is.na(b)], "\n")
      else
        cat("Unknown.\n")
    }
  }
}
\keyword{graphs}
\keyword{models}
\keyword{multivariate}

