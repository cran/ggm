\name{ggm}
\alias{ggm}
\title{
The package \code{ggm}: summary information
}
\description{
This package provides functions for defining, manipulating and fitting graphical Markov models with mixed graphs. It is intended as a contribution to the gR-project described by Lauritzen (2002).

For a tutorial illustrating the new functions in the package `ggm` that deal with ancestral, summary and ribbonless graphs see Sadeghi and Marchetti (2012) in the references.
}
\section{Functions}{
The main functions can be classified as follows.
\itemize{
\item Functions for defining graphs (undirected, directed acyclic,
ancestral and summary graphs): \code{\link{UG}}, \code{\link{DAG}}, \code{\link{makeMG}}, \code{\link{grMAT}};
\item Functions for doing graph operations (parents, boundary,
cliques, connected components, fundamental cycles, d-separation, m-separation):
\code{\link{pa}}, \code{\link{bd}}, \code{\link{cliques}},
\code{\link{conComp}}, \code{\link{fundCycles}};
\item Functions for testing independence statements and generating maximal graphs
from non-maximal graphs:
\code{\link{dSep}}, \code{\link{msep}}, \code{\link{Max}};
\item Function for finding covariance and concentration graphs
induced by marginalization and conditioning:
\code{\link{inducedCovGraph}}, \code{\link{inducedConGraph}};
\item Functions for finding multivariate regression graphs and chain
graphs induced by marginalization and conditioning:
\code{\link{inducedRegGraph}},  \code{\link{inducedChainGraph}},
\code{\link{inducedDAG}};
\item Functions for finding stable mixed graphs (ancestral, summary and
ribbonless) after marginalization and conditioning:
\code{\link{AG}},  \code{\link{SG}}, \code{\link{RG}};
\item Functions for fitting by ML Gaussian DAGs, concentration graphs,
covariance graphs and ancestral graphs:
\code{\link{fitDag}}, \code{\link{fitConGraph}},
\code{\link{fitCovGraph}}, \code{\link{fitAncestralGraph}};
\item Functions for testing several conditional independences:\code{\link{shipley.test}};
\item Functions for checking global identification
of DAG Gaussian  models with one latent variable (Stanghellini-Vicard's condition for concentration graphs, new sufficient conditions for DAGs):
\code{\link{isGident}}, \code{\link{checkIdent}};
\item Functions for fitting Gaussian DAG models with one latent variable:
\code{\link{fitDagLatent}};
\item Functions for testing Markov equivalences and generating Markov
equivalent graphs of specific types:
\code{\link{MarkEqRcg}},  \code{\link{MarkEqMag}},
\code{\link{RepMarDAG}}, \code{\link{RepMarUG}}, \code{\link{RepMarBG}}.
}
}
\section{Authors}{
Giovanni M. Marchetti, Dipartimento di Statistica, Informatica, Applicazioni 'G. Parenti'. University of Florence, Italy

Mathias Drton,  Department of Statistics, University of Washington, USA

Kayvan Sadeghi, Department of Statistics, Carnegie Mellon University, USA
}
\section{Acknowledgements}{
Many thanks to Fulvia Pennoni for testing some of
the functions, to Elena Stanghellini for discussion and
examples and to Claus Dethlefsen and Jens Henrik Badsberg for
suggestions and corrections. The function \code{fitConGraph} was corrected by
Ilaria Carobbi. Helpful discussions with Steffen Lauritzen and Nanny Wermuth,
are gratefully acknowledged. Thanks also to Michael Perlman,
Thomas Richardson and David Edwards.

Giovanni Marchetti has been supported by MIUR, Italy, under grant
scheme PRIN 2002, and Mathias Drton has been supported by NSF grant
DMS-9972008 and University of Washington RRF grant 65-3010.
}
\references{
Lauritzen, S. L. (2002). gRaphical Models in R.
\emph{R News}, 3(2)39.

Sadeghi, K. and Marchetti, G.M. (2012). Graphical Markov models with mixed graphs in R. \emph{The R Journal},  4(2):65-73. \url{https://journal.r-project.org/archive/2012/RJ-2012-015/RJ-2012-015.pdf}
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

