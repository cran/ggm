\name{ggm}
\alias{ggm}
\title{
The package `ggm': summary information
}
\description{
This package provides functions for defining and manipulating graphs
and for fitting graphical Gaussian models.
}
\section{Functions}{
The main functions can be classified as follows.
\itemize{
\item Functions for defining graphs (undirected, directed acyclic, 
ancestral graphs): \code{\link{UG}}, \code{\link{DAG}}, \code{\link{makeAG}};
\item Functions for doing graph operations (parents, boundary, 
cliques, connected components, fundamental cycles, d-separation):
\code{\link{pa}}, \code{\link{bd}}, \code{\link{cliques}},
\code{\link{conComp}}, \code{\link{fundCycles}}, \code{\link{dSep}};
\item Function for finding covariance and concentration graphs
induced by marginalization and conditioning:
\code{\link{inducedCovGraph}}, \code{\link{inducedConGraph}};
\item Functions for fitting by ML Gaussian DAGs, concentration graphs, 
covariance graphs and ancestral graphs:
\code{\link{fitDag}}, \code{\link{fitConGraph}},
\code{\link{fitCovGraph}}, \code{\link{fitAncestralGraph}};
\item Functions for testing several conditional independences:
\code{\link{shipley.test}};
\item Functions for checking global identification 
of DAG Gaussian  models with one 
latent variable (Stanghellini-Vicard's condition for concentration graphs, 
new sufficient conditions for DAGs):
\code{\link{isGident}}, \code{\link{checkIdent}};
\item Functions for fitting Gaussian DAG models with one latent variable:
\code{\link{fitDagLatent}}.
}
The package is intended as a contribution to the gR-project 
derscribed by Lauritzen (2002).
}
\section{Authors}{
Giovanni M. Marchetti, Dipartimento di Statistica ``G. Parenti''.
Universit\`a di Firenze, Italy; 
 
Mathias Drton,  Department of Statistics, University of Washington, USA.
}
\section{Acknowledgements}{
Many thanks to Fulvia Pennoni for testing some of 
the functions, and to Elena Stanghellini for discussion and
examples. Thanks also to Michael Perlman and Thomas Richardson.

Giovanni Marchetti has been supported by MIUR, Italy, under grant 
scheme PRIN 2002, and Mathias Drton has been supported by NSF grant 
DMS-9972008 and University of Washington RRF grant 65-3010.
}
\references{
Lauritzen, S. L. (2002). gRaphical Models in R.
\emph{R News}, 3(2)39.
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
 
