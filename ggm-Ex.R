pkgname <- "ggm"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ggm')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("AG")
### * AG

flush(stderr()); flush(stdout())

### Name: AG
### Title: Ancestral graph
### Aliases: AG
### Keywords: graphs ancestral graph directed acyclic graph marginalization
###   and conditioning

### ** Examples
                
##The adjacency matrix of a DAG    
ex<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
             0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0),16,16,byrow=TRUE)
M <- c(3,5,6,15,16)
C <- c(4,7)
AG(ex, M, C, plot = TRUE)



cleanEx()
nameEx("DAG")
### * DAG

flush(stderr()); flush(stdout())

### Name: DAG
### Title: Directed acyclic graphs (DAGs)
### Aliases: DAG
### Keywords: graphs models multivariate

### ** Examples

## A Markov chain
DAG(y ~ x, x ~ z, z ~ u)

## Another DAG
DAG(y ~ x + z + u, x ~ u, z ~ u)

## A DAG with an isolated node
DAG(v ~ v, y ~ x + z, z ~ w + u)

## There can be repetitions
DAG(y ~ x + u + v, y ~ z, u ~ v + z)

## Interactions are ignored
DAG(y ~ x*z + z*v, x ~ z)

## A cyclic graph returns an error!
## Not run: DAG(y ~ x, x ~ z, z ~ y)

## The order can be changed
DAG(y ~ z, y ~ x + u + v,  u ~ v + z)

## If you want to order the nodes (topological sort of the DAG)
DAG(y ~ z, y ~ x + u + v,  u ~ v + z, order=TRUE)



cleanEx()
nameEx("DG")
### * DG

flush(stderr()); flush(stdout())

### Name: DG
### Title: Directed graphs
### Aliases: DG
### Keywords: graphs directed graph models multivariate

### ** Examples

## A DAG
DG(y ~ x, x ~ z, z ~ u)

## A cyclic directed graph
DG(y ~ x, x ~ z, z ~ y)

## A graph with two arrows between two nodes
DG(y ~ x, x ~ y)

## There can be isolated nodes
DG(y ~ x, x ~ x)



cleanEx()
nameEx("In")
### * In

flush(stderr()); flush(stdout())

### Name: In
### Title: Indicator matrix
### Aliases: In
### Keywords: array algebra graphs multivariate

### ** Examples

## A simple way to find the overall induced concentration graph
## The DAG on p. 198 of Cox & Wermuth (1996)
amat <- DAG(y1 ~ y2 + y3, y3 ~ y5, y4 ~ y5)
A <- edgematrix(amat)
In(crossprod(A))



cleanEx()
nameEx("InducedGraphs")
### * InducedGraphs

flush(stderr()); flush(stdout())

### Name: InducedGraphs
### Title: Graphs induced by marginalization or conditioning
### Aliases: inducedCovGraph inducedConGraph inducedRegGraph
###   inducedChainGraph inducedDAG InducedGraphs
### Keywords: graphs models multivariate

### ** Examples

## Define a DAG
dag <- DAG(a ~ x, c ~ b+d, d~ x)
dag
## Induced covariance graph of a, b, d given the empty set.
inducedCovGraph(dag, sel=c("a", "b", "d"), cond=NULL)

## Induced concentration graph of a, b, c given x
inducedConGraph(dag, sel=c("a", "b", "c"), cond="x")

## Overall covariance graph
inducedCovGraph(dag)

## Overall concentration graph
inducedConGraph(dag)

## Induced covariance graph of x, b, d given c, x.
inducedCovGraph(dag, sel=c("a", "b", "d"), cond=c("c", "x"))

## Induced concentration graph of a, x, c given d, b.
inducedConGraph(dag, sel=c("a", "x", "c"), cond=c("d", "b"))

## The DAG on p. 198 of Cox & Wermuth (1996)
dag <- DAG(y1~ y2 + y3, y3 ~ y5, y4 ~ y5)

## Cf. figure 8.7 p. 203 in Cox & Wermuth (1996)
inducedCovGraph(dag, sel=c("y2", "y3", "y4", "y5"), cond="y1")
inducedCovGraph(dag, sel=c("y1", "y2", "y4", "y5"), cond="y3")
inducedCovGraph(dag, sel=c("y1", "y2", "y3", "y4"), cond="y5")

## Cf. figure 8.8 p. 203 in Cox & Wermuth (1996)
inducedConGraph(dag, sel=c("y2", "y3", "y4", "y5"), cond="y1")
inducedConGraph(dag, sel=c("y1", "y2", "y4", "y5"), cond="y3")
inducedConGraph(dag, sel=c("y1", "y2", "y3", "y4"), cond="y5")

## Cf. figure 8.9 p. 204 in Cox & Wermuth (1996)
inducedCovGraph(dag, sel=c("y2", "y3", "y4", "y5"), cond=NULL)
inducedCovGraph(dag, sel=c("y1", "y2", "y4", "y5"), cond=NULL)
inducedCovGraph(dag, sel=c("y1", "y2", "y3", "y4"), cond=NULL)

## Cf. figure 8.10 p. 204 in Cox & Wermuth (1996)
inducedConGraph(dag, sel=c("y2", "y3", "y4", "y5"), cond=NULL)
inducedConGraph(dag, sel=c("y1", "y2", "y4", "y5"), cond=NULL)
inducedConGraph(dag, sel=c("y1", "y2", "y3", "y4"), cond=NULL)

## An induced regression graph
dag2 = DAG(Y ~ X+U, W ~ Z+U)
inducedRegGraph(dag2, sel="W",  cond=c("Y", "X", "Z"))

## An induced DAG
inducedDAG(dag2, order=c("X","Y","Z","W"))

## An induced multivariate regression graph
inducedRegGraph(dag2, sel=c("Y", "W"), cond=c("X", "Z"))

## An induced chain graph with LWF interpretation
dag3 = DAG(X~W, W~Y, U~Y+Z)
cc = list(c("W", "U"), c("X", "Y", "Z"))
inducedChainGraph(dag3, cc=cc, type="LWF")

## ... with AMP interpretation
inducedChainGraph(dag3, cc=cc, type="AMP")

## ... with multivariate regression interpretation
cc= list(c("U"), c("Z", "Y"), c("X", "W"))
inducedChainGraph(dag3, cc=cc, type="MRG")



cleanEx()
nameEx("MAG")
### * MAG

flush(stderr()); flush(stdout())

### Name: MAG
### Title: Maximal ancestral graph
### Aliases: MAG
### Keywords: ancestral graph directed acyclic graph marginalization and
###   conditioning maximality of graphs

### ** Examples

ex<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ##The adjacency matrix of a DAG
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
             0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0), 16, 16, byrow = TRUE)
M <- c(3,5,6,15,16)
C <- c(4,7)
MAG(ex, M, C, plot=TRUE)
###################################################
H <- matrix(c(0,100,1,0,100,0,100,0,0,100,0,100,0,1,100,0),4,4)
Max(H)



cleanEx()
nameEx("MRG")
### * MRG

flush(stderr()); flush(stdout())

### Name: MRG
### Title: Maximal ribbonless graph
### Aliases: MRG
### Keywords: graphs directed acyclic graph marginalisation and
###   conditioning maximality of graphs MC graph ribbonless graph

### ** Examples

ex <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ##The adjacency matrix of a DAG
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
               0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
               0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,
               1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
               0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0),16,16, byrow = TRUE)
M <- c(3,5,6,15,16)
C <- c(4,7)
MRG(ex, M, C, plot = TRUE)
###################################################
H <- matrix(c( 0, 100,   1,   0,
	         100,   0, 100,   0,
	           0, 100,   0, 100,
	           0,   1, 100,   0), 4,4)
Max(H)



cleanEx()
nameEx("MSG")
### * MSG

flush(stderr()); flush(stdout())

### Name: MSG
### Title: Maximal summary graph
### Aliases: MSG
### Keywords: graphs directed acyclic graph marginalisation and
###   conditioning maximality of graphs summary graph

### ** Examples

ex<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ##The adjacency matrix of a DAG
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
             0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0), 16, 16, byrow=TRUE)
M <- c(3,5,6,15,16)
C <- c(4,7)
MSG(ex,M,C,plot=TRUE)
###################################################
H<-matrix(c(0,100,1,0,100,0,100,0,0,100,0,100,0,1,100,0),4,4)
Max(H)



cleanEx()
nameEx("MarkEqMag")
### * MarkEqMag

flush(stderr()); flush(stdout())

### Name: MarkEqMag
### Title: Markov equivalence of maximal ancestral graphs
### Aliases: MarkEqMag
### Keywords: graphs Markov equivalence maximal ancestral graphs
###   multivariate

### ** Examples

H1<-matrix(  c(0,100,  0,  0,
	         100,  0,100,  0,
               0,100,  0,100,
               0,  1,100,  0), 4, 4)
H2<-matrix(c(0,0,0,0,1,0,100,0,0,100,0,100,0,1,100,0),4,4)
H3<-matrix(c(0,0,0,0,1,0,0,0,0,1,0,100,0,1,100,0),4,4)
MarkEqMag(H1,H2)
MarkEqMag(H1,H3)
MarkEqMag(H2,H3)



cleanEx()
nameEx("MarkEqRcg")
### * MarkEqRcg

flush(stderr()); flush(stdout())

### Name: MarkEqRcg
### Title: Markov equivalence for regression chain graphs.
### Aliases: MarkEqRcg
### Keywords: graphs bidirected graph directed acyclic graph Markov
###   equivalence regression chain graph undirected graph multivariate

### ** Examples

H1<-matrix(c(0,100,0,0,0,100,0,100,0,0,0,100,0,0,0,1,0,0,0,100,0,0,1,100,0),5,5)
H2<-matrix(c(0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,100,0,0,1,100,0),5,5)
H3<-matrix(c(0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0),5,5)
#MarkEqRcg(H1,H2)
#MarkEqRcg(H1,H3)
#MarkEqRcg(H2,H3)



cleanEx()
nameEx("Max")
### * Max

flush(stderr()); flush(stdout())

### Name: Max
### Title: Maximisation for graphs
### Aliases: Max
### Keywords: graphs loopless mixed graph m-separation maximality

### ** Examples

H <- matrix(c(  0,100,  1,  0,
	          100,  0,100,  0,
	            0,100,  0,100,
	            0,  1,100,  0), 4, 4)
Max(H)



cleanEx()
nameEx("RG")
### * RG

flush(stderr()); flush(stdout())

### Name: RG
### Title: Ribbonless graph
### Aliases: RG
### Keywords: graphs directed acyclic graph marginalisation and
###   conditioning MC graph ribbonless graph

### ** Examples

	ex <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ##The adjacency matrix of a DAG
	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	               1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	               0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
	               0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	               0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
	               0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
	               0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,
	               1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	               0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,
	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	               1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
	               0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0),16,16, byrow = TRUE)

M<-c(3,5,6,15,16)
C<-c(4,7)
RG(ex,M,C,plot=TRUE)



cleanEx()
nameEx("RepMarBG")
### * RepMarBG

flush(stderr()); flush(stdout())

### Name: RepMarBG
### Title: Representational Markov equivalence to bidirected graphs.
### Aliases: RepMarBG
### Keywords: graphs bidirected graph Markov equivalence maximal ancestral
###   graph representational Markov equivalence

### ** Examples

H<-matrix(c(0,10,0,0,10,0,0,0,0,1,0,100,0,0,100,0),4,4)
RepMarBG(H)



cleanEx()
nameEx("RepMarDAG")
### * RepMarDAG

flush(stderr()); flush(stdout())

### Name: RepMarDAG
### Title: Representational Markov equivalence to directed acyclic graphs.
### Aliases: RepMarDAG
### Keywords: graphs bidirected graph Markov equivalence maximal ancestral
###   graph representational Markov equivalence

### ** Examples

H<-matrix(c(0,10,0,0,10,0,0,0,0,1,0,100,0,0,100,0),4,4)
RepMarBG(H)



cleanEx()
nameEx("RepMarUG")
### * RepMarUG

flush(stderr()); flush(stdout())

### Name: RepMarUG
### Title: Representational Markov equivalence to undirected graphs.
### Aliases: RepMarUG
### Keywords: graphs bidirected graph Markov equivalence maximal ancestral
###   graph representational Markov equivalence

### ** Examples

H<-matrix(c(0,10,0,0,10,0,0,0,0,1,0,100,0,0,100,0),4,4)
RepMarUG(H)



cleanEx()
nameEx("SG")
### * SG

flush(stderr()); flush(stdout())

### Name: SG
### Title: summary graph
### Aliases: SG
### Keywords: graphs directed acyclic graph marginalization and
###   conditioning summary graph

### ** Examples

	ex <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ##The adjacency matrix of a DAG
	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	               1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	               0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
	               0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	               0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
	               0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
	               0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,
	               1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	               0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,
	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	               1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
	               0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0),16,16, byrow = TRUE)
M <- c(3,5,6,15,16)
C <- c(4,7)
SG(ex, M, C, plot = TRUE)



cleanEx()
nameEx("SimpleGraphOperations")
### * SimpleGraphOperations

flush(stderr()); flush(stdout())

### Name: Simple Graph Operations
### Title: Simple graph operations
### Aliases: bd ch pa
### Keywords: graphs models multivariate

### ** Examples

## find boundary of a subset of nodes of a DAG
G <- DAG(y ~ x+b+a, b~a, x~a)
bd("b", G)
bd(c("b", "x"), G)
bd("x", G)
bd(c("x","b"), G)
## find boundary of a subset of nodes of an UG
G <- UG(~ y*x*z + z*h*v)
bd("z", G)
bd(c("y", "x"), G)
bd("v", G)
bd(c("x","v"), G)
## children of a subset of nodes of a DAG
G <- DAG(y ~ x+b+a, b~a, x~a)
ch("b", G)
ch(c("b", "x"), G)
ch("x", G)
ch(c("a","x"), G)
## parents of a subset of nodes of a DAG
pa("b", G)
pa(c("b", "x"), G)
pa("x", G)
pa(c("x","b"), G)



cleanEx()
nameEx("UG")
### * UG

flush(stderr()); flush(stdout())

### Name: UG
### Title: Defining an undirected graph (UG)
### Aliases: UG
### Keywords: graphs models multivariate

### ** Examples

## X independent of Y given Z
UG(~ X*Z + Y*Z)

# The saturated model
UG(~ X*Y*Z)

## The model without three-way interactions has the same graph
UG(~ X*Y + Y*Z + Z*X)
UG(~ (X + Y + Z)^2)

## Butterfly model defined from the cliques
UG(~ mec*vec*alg + alg*ana*sta)

## Some isolated nodes
UG(~x*y*z + a + b) 



cleanEx()
nameEx("adjMatrix")
### * adjMatrix

flush(stderr()); flush(stdout())

### Name: adjMatrix
### Title: Adjacency matrix of a graph
### Aliases: adjMatrix
### Keywords: array algebra graphs multivariate

### ** Examples

amat <- DAG(y ~ x+z, z~u+v)
E <- edgematrix(amat)
adjMatrix(E)



cleanEx()
nameEx("allEdges")
### * allEdges

flush(stderr()); flush(stdout())

### Name: allEdges
### Title: All edges of a graph
### Aliases: allEdges
### Keywords: graphs models multivariate

### ** Examples

## A UG graph
allEdges(UG(~ y*v*k +v*k*d+y*d))

## A DAG
allEdges(DAG(u~h+o+p, h~o, o~p))



cleanEx()
nameEx("anger")
### * anger

flush(stderr()); flush(stdout())

### Name: anger
### Title: Anger data
### Aliases: anger
### Keywords: datasets

### ** Examples
 
# Fit a chordless 4-cycle model 
data(anger) 
G = UG(~ Y*X + X*Z + Z*U + U*Y)
fitConGraph(G,anger, 684) 



cleanEx()
nameEx("basiSet")
### * basiSet

flush(stderr()); flush(stdout())

### Name: basiSet
### Title: Basis set of a DAG
### Aliases: basiSet
### Keywords: graphs models multivariate

### ** Examples

## See Shipley (2000), Figure 2, p. 213
A <- DAG(x5~ x3+x4, x3~ x2, x4~x2, x2~ x1)
basiSet(A)



cleanEx()
nameEx("bfsearch")
### * bfsearch

flush(stderr()); flush(stdout())

### Name: bfsearch
### Title: Breadth first search
### Aliases: bfsearch
### Keywords: graphs models multivariate

### ** Examples

## Finding a spanning tree of the butterfly graph
bfsearch(UG(~ a*b*o + o*u*j))
## Starting from another node
bfsearch(UG(~ a*b*o + o*u*j), v=3)



cleanEx()
nameEx("blkdiag")
### * blkdiag

flush(stderr()); flush(stdout())

### Name: blkdiag
### Title: Block diagonal matrix
### Aliases: blkdiag
### Keywords: matrix

### ** Examples

X <- c(1,1,2,2); Z <- c(10, 20, 30, 40); A <- factor(c(1,2,2,2))
blkdiag(model.matrix(~X+Z), model.matrix(~A))



cleanEx()
nameEx("blodiag")
### * blodiag

flush(stderr()); flush(stdout())

### Name: blodiag
### Title: Block diagonal matrix
### Aliases: blodiag
### Keywords: matrix

### ** Examples

blodiag(1:10, blo = c(2, 3, 5)) 
blodiag(1:10, blo = c(3,4,0,1))



cleanEx()
nameEx("checkIdent")
### * checkIdent

flush(stderr()); flush(stdout())

### Name: checkIdent
### Title: Identifiability of a model with one latent variable
### Aliases: checkIdent
### Keywords: graphs models multivariate

### ** Examples

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
      if(TRUE %in% b)
        cat("Node", i, names(b)[!is.na(b)], "\n")
      else
        cat("Unknown.\n")
    }
  }



cleanEx()
nameEx("cmpGraph")
### * cmpGraph

flush(stderr()); flush(stdout())

### Name: cmpGraph
### Title: The complementary graph
### Aliases: cmpGraph
### Keywords: graphs models multivariate

### ** Examples

## A chordless four-cycle
four <- UG(~ a*b + b*d + d*e + e*a)
four
cmpGraph(four)



cleanEx()
nameEx("conComp")
### * conComp

flush(stderr()); flush(stdout())

### Name: conComp
### Title: Connectivity components
### Aliases: conComp
### Keywords: graphs models multivariate

### ** Examples

## three connected components
conComp(UG(~a*c+c*d+e+g*o*u))
## a connected graph
conComp(UG(~ a*b+b*c+c*d+d*a))



cleanEx()
nameEx("correlations")
### * correlations

flush(stderr()); flush(stdout())

### Name: correlations
### Title: Marginal and partial correlations
### Aliases: correlations
### Keywords: array graphs models multivariate

### ** Examples

## See Table 6.1 in Cox & Wermuth (1996)
data(glucose)
correlations(glucose)



cleanEx()
nameEx("cycleMatrix")
### * cycleMatrix

flush(stderr()); flush(stdout())

### Name: cycleMatrix
### Title: Fundamental cycles
### Aliases: cycleMatrix
### Keywords: graphs models multivariate

### ** Examples

## Three cycles
cycleMatrix(UG(~a*b*d+d*e+e*a*f))
## No cycle
 cycleMatrix(UG(~a*b))
## two cycles: the first is even and the second is odd
cm <- cycleMatrix(UG(~a*b+b*c+c*d+d*a+a*u*v))
apply(cm, 1, sum)



cleanEx()
nameEx("dSep")
### * dSep

flush(stderr()); flush(stdout())

### Name: dSep
### Title: d-separation
### Aliases: dSep
### Keywords: graphs models multivariate

### ** Examples

## Conditioning on a transition node
dSep(DAG(y ~ x, x ~ z), first="y", second="z", cond = "x")
## Conditioning on a collision node (collider)
dSep(DAG(y ~ x, y ~ z), first="x", second="z", cond = "y")
## Conditioning on a source node
dSep(DAG(y ~ x, z ~ x), first="y", second="z", cond = "x")
## Marginal independence
dSep(DAG(y ~ x, y ~ z), first="x", second="z", cond = NULL)
## The DAG defined on p.~47 of Lauritzen (1996)
dag <- DAG(g ~ x, h ~ x+f, f ~ b, x ~ l+d, d ~ c, c ~ a, l ~ y, y ~ b)
dSep(dag, first="a", second="b", cond=c("x", "y"))
dSep(dag, first="a", second=c("b", "d"), cond=c("x", "y"))



cleanEx()
nameEx("derived")
### * derived

flush(stderr()); flush(stdout())

### Name: derived
### Title: Data on blood pressure body mass and age
### Aliases: derived
### Keywords: datasets

### ** Examples

# A DAG model with a latent variable U
G = DAG(Y ~ Z + U, X ~ U + W, Z ~ W)

data(derived)

# The model fitted using the derived variables
out = fitDagLatent(G, derived$S, n = 44, latent = "U")

# An ancestral graph model marginalizing over U
H = AG(G, M = "U")

# The ancestral graph model fitted obtaining the 
# same result
out2 = fitAncestralGraph(H, derived$S, n = 44)



cleanEx()
nameEx("diagv")
### * diagv

flush(stderr()); flush(stdout())

### Name: diagv
### Title: Matrix product with a diagonal matrix
### Aliases: diagv
### Keywords: matrix

### ** Examples

v <- 1:1000
M <- matrix(runif(3000), 1000, 3)
dim(diagv(v, M))



cleanEx()
nameEx("drawGraph")
### * drawGraph

flush(stderr()); flush(stdout())

### Name: drawGraph
### Title: Drawing a graph with a simple point and click interface.
### Aliases: drawGraph
### Keywords: graphs hplot iplot

### ** Examples

## A directed acyclic graph
d <- DAG(y1 ~ y2+y6, y2 ~ y3, y3 ~ y5+y6, y4 ~ y5+y6)
## Not run: drawGraph(d)

## An undirected graph
g <- UG(~giova*anto*armo + anto*arj*sara) 
## Not run: drawGraph(d)

## An ancestral graph
ag <- makeMG(ug=UG(~y0*y1), dg=DAG(y4~y2, y2~y1), bg=UG(~y2*y3+y3*y4))
drawGraph(ag, adjust = FALSE)
drawGraph(ag, adjust = FALSE)

## A more complex example with coordinates: the UNIX evolution
xy <-
structure(c(5, 15, 23, 25, 26, 17, 8, 6, 6, 7, 39, 33, 23, 49, 
19, 34, 13, 29, 50, 68, 70, 86, 89, 64, 81, 45, 64, 49, 64, 87, 
65, 65, 44, 37, 64, 68, 73, 85, 83, 95, 84, 0, 7, 15, 27, 44, 
37, 36, 20, 51, 65, 44, 64, 59, 73, 69, 78, 81, 90, 97, 89, 72, 
85, 74, 62, 68, 59, 52, 48, 43, 50, 34, 21, 18, 5, 1, 10, 2, 
11, 2, 1, 44), .Dim = c(41, 2), .Dimnames = list(NULL, c("x", 
"y")))
Unix <- DAG(
                SystemV.3 ~ SystemV.2,
                SystemV.2 ~ SystemV.0,
                SystemV.0 ~ TS4.0,
                TS4.0 ~ Unix.TS3.0 + Unix.TS.PP + CB.Unix.3,
                PDP11.SysV ~ CB.Unix.3,
                CB.Unix.3 ~ CB.Unix.2,
                CB.Unix.2 ~ CB.Unix.1,
                Unix.TS.PP ~ CB.Unix.3,
                Unix.TS3.0 ~ Unix.TS1.0 + PWB2.0 + USG3.0 + Interdata,
                USG3.0 ~ USG2.0,
                PWB2.0 ~ Interdata + PWB1.2,
                USG2.0 ~ USG1.0,
                CB.Unix.1 ~ USG1.0,
                PWB1.2 ~ PWB1.0,
                USG1.0 ~ PWB1.0,
                PWB1.0 ~ FifthEd,
                SixthEd ~ FifthEd,
                LSX ~ SixthEd,
                MiniUnix ~ SixthEd,
                Interdata ~ SixthEd,
                Wollongong ~ SixthEd,
                SeventhEd ~ Interdata,
                BSD1 ~ SixthEd,
                Xenix ~ SeventhEd,
                V32 ~ SeventhEd,
                Uniplus ~ SeventhEd,
                BSD3 ~ V32,
                BSD2 ~ BSD1,
                BSD4 ~ BSD3,
                BSD4.1 ~ BSD4,
                EigthEd ~ SeventhEd + BSD4.1,
                NinethEd ~ EigthEd,
                Ultrix32 ~ BSD4.2,
                BSD4.2 ~ BSD4.1,
                BSD4.3 ~ BSD4.2,
                BSD2.8 ~ BSD4.1 + BSD2,
                BSD2.9 ~ BSD2.8,
                Ultrix11 ~ BSD2.8 + V7M + SeventhEd,
                V7M ~ SeventhEd
                )
drawGraph(Unix, coor=xy, adjust=FALSE)
# dev.print(file="unix.fig", device=xfig) # Edit the graph with Xfig



cleanEx()
nameEx("edgematrix")
### * edgematrix

flush(stderr()); flush(stdout())

### Name: edgematrix
### Title: Edge matrix of a graph
### Aliases: edgematrix
### Keywords: array algebra graphs multivariate

### ** Examples

amat <- DAG(y ~ x+z, z~u+v)
amat
edgematrix(amat)
edgematrix(amat, inv=TRUE)



cleanEx()
nameEx("essentialGraph")
### * essentialGraph

flush(stderr()); flush(stdout())

### Name: essentialGraph
### Title: Essential graph
### Aliases: essentialGraph
### Keywords: graphs models multivariate

### ** Examples

dag = DAG(U ~ Y+Z, Y~X, Z~X)
essentialGraph(dag)



cleanEx()
nameEx("findPath")
### * findPath

flush(stderr()); flush(stdout())

### Name: findPath
### Title: Finding paths
### Aliases: findPath
### Keywords: graphs

### ** Examples

## A (single) path on a spanning tree
findPath(bfsearch(UG(~ a*b*c + b*d + d*e+ e*c))$tree, st=1, en=5)



cleanEx()
nameEx("fitAncestralGraph")
### * fitAncestralGraph

flush(stderr()); flush(stdout())

### Name: fitAncestralGraph
### Title: Fitting of Gaussian Ancestral Graph Models
### Aliases: fitAncestralGraph
### Keywords: graphs models ancestral graph multivariate

### ** Examples

## A covariance matrix
"S" <- structure(c(2.93, -1.7, 0.76, -0.06,
                  -1.7, 1.64, -0.78, 0.1,
                   0.76, -0.78, 1.66, -0.78,
                  -0.06, 0.1, -0.78, 0.81), .Dim = c(4,4),
                 .Dimnames = list(c("y", "x", "z", "u"), c("y", "x", "z", "u")))
## The following should give the same fit.   
## Fit an ancestral graph y -> x <-> z <- u
fitAncestralGraph(ag1 <- makeMG(dg=DAG(x~y,z~u), bg = UG(~x*z)), S, n=100)

## Fit an ancestral graph y <-> x <-> z <-> u
fitAncestralGraph(ag2 <- makeMG(bg= UG(~y*x+x*z+z*u)), S, n=100)

## Fit the same graph with fitCovGraph
fitCovGraph(ag2, S, n=100)    

## Another example for the mathematics marks data

data(marks)
S <- var(marks)
mag1 <- makeMG(bg=UG(~mechanics*vectors*algebra+algebra*analysis*statistics))
fitAncestralGraph(mag1, S, n=88)

mag2 <- makeMG(ug=UG(~mechanics*vectors+analysis*statistics),
               dg=DAG(algebra~mechanics+vectors+analysis+statistics))
fitAncestralGraph(mag2, S, n=88) # Same fit as above



cleanEx()
nameEx("fitConGraph")
### * fitConGraph

flush(stderr()); flush(stdout())

### Name: fitConGraph
### Title: Fitting a Gaussian concentration graph model
### Aliases: fitConGraph
### Keywords: graphs models multivariate

### ** Examples

## A model for the mathematics marks (Whittaker, 1990)
data(marks)
## A butterfly concentration graph  
G <- UG(~ mechanics*vectors*algebra + algebra*analysis*statistics)
fitConGraph(G, cov(marks), nrow(marks))   
## Using the cliques

cl = list(c("mechanics", "vectors",   "algebra"), c("algebra", "analysis" ,  "statistics")) 
fitConGraph(G, S = cov(marks), n = nrow(marks), cli = cl) 



cleanEx()
nameEx("fitCovGraph")
### * fitCovGraph

flush(stderr()); flush(stdout())

### Name: fitCovGraph
### Title: Fitting of Gaussian covariance graph models
### Aliases: fitCovGraph
### Keywords: graphs models multivariate

### ** Examples

## Correlations among four strategies to cope with stress for 
## 72 students. Cox & Wermuth (1996), p. 73.

data(stress)

## A chordless 4-cycle covariance graph
G <- UG(~ Y*X + X*U + U*V + V*Y)

fitCovGraph(G, S = stress, n=72)
fitCovGraph(G, S = stress, n=72, alg="dual")



cleanEx()
nameEx("fitDag")
### * fitDag

flush(stderr()); flush(stdout())

### Name: fitDag
### Title: Fitting of Gaussian DAG models
### Aliases: fitDag
### Keywords: graphs models multivariate

### ** Examples

dag <- DAG(y ~ x+u, x ~ z, z ~ u)
"S" <- structure(c(2.93, -1.7, 0.76, -0.06,
                   -1.7, 1.64, -0.78, 0.1,
                    0.76, -0.78, 1.66, -0.78,
                    -0.06, 0.1, -0.78, 0.81), .Dim = c(4,4),
         .Dimnames = list(c("y", "x", "z", "u"), c("y", "x", "z", "u")))
fitDag(dag, S, 200)



cleanEx()
nameEx("fitDagLatent")
### * fitDagLatent

flush(stderr()); flush(stdout())

### Name: fitDagLatent
### Title: Fitting Gaussian DAG models with one latent variable
### Aliases: fitDagLatent
### Keywords: graphs models multivariate

### ** Examples

## data from Joreskog and Goldberger (1975)
V <- matrix(c(1,     0.36,   0.21,  0.10,  0.156, 0.158,
              0.36,  1,      0.265, 0.284, 0.192, 0.324,
              0.210, 0.265,  1,     0.176, 0.136, 0.226,
              0.1,   0.284,  0.176, 1,     0.304, 0.305, 
              0.156, 0.192,  0.136, 0.304, 1,     0.344,
              0.158, 0.324,  0.226, 0.305, 0.344, 1),     6,6)
nod <- c("y1", "y2", "y3", "x1", "x2", "x3")
dimnames(V) <- list(nod,nod)
dag <- DAG(y1 ~ z, y2 ~ z, y3 ~ z, z ~ x1 + x2 + x3, x1~x2+x3, x2~x3) 
fitDagLatent(dag, V, n=530, latent="z", seed=4564)
fitDagLatent(dag, V, n=530, latent="z", norm=2, seed=145)



cleanEx()
nameEx("fitmlogit")
### * fitmlogit

flush(stderr()); flush(stdout())

### Name: fitmlogit
### Title: Multivariate logistic models
### Aliases: fitmlogit
### Keywords: multivariate logistic model

### ** Examples
    
data(surdata)                     
out1 <- fitmlogit(A ~X, B ~ Z, cbind(A, B) ~ X*Z, data = surdata)     
out1$beta
out2 <- fitmlogit(A ~X, B ~ Z, cbind(A, B) ~ 1, data = surdata)        
out2$beta



cleanEx()
nameEx("fundCycles")
### * fundCycles

flush(stderr()); flush(stdout())

### Name: fundCycles
### Title: Fundamental cycles
### Aliases: fundCycles
### Keywords: graphs models multivariate

### ** Examples

## Three fundamental cycles
fundCycles(UG(~a*b*d + d*e + e*a*f))



cleanEx()
nameEx("glucose")
### * glucose

flush(stderr()); flush(stdout())

### Name: glucose
### Title: Glucose control
### Aliases: glucose
### Keywords: datasets

### ** Examples

data(glucose)
## See Cox & Wermuth (1996), Figure 6.3 p. 140
coplot(Y ~ W | A, data=glucose)



cleanEx()
nameEx("grMAT")
### * grMAT

flush(stderr()); flush(stdout())

### Name: grMAT
### Title: Graph to adjacency matrix
### Aliases: grMAT
### Keywords: graphs adjacency matrix mixed graph vector

### ** Examples

## Generating the adjacency matrix from an igraph object
exdag <- graph.formula(v7+-v8-+v5+-v6, v5-+v1+-v4-+v3, v1+-v2+-v3)
grMAT(exdag)

## Generating the adjacency matrix from a vector
exvec <-c ('b',1,2,'b',1,14,'a',9,8,'l',9,11,'a',10,8,
           'a',11,2,'a',11,10,'a',12,1,'b',12,14,'a',13,10,'a',13,12)
grMAT(exvec)



cleanEx()
nameEx("isADMG")
### * isADMG

flush(stderr()); flush(stdout())

### Name: isADMG
### Title: Acyclic directed mixed graphs
### Aliases: isADMG
### Keywords: graphs ancestral graph mixed graph models multivariate

### ** Examples

	## Examples from Richardson and Spirtes (2002)
	a1 <- makeMG(dg=DAG(a~b, b~d, d~c), bg=UG(~a*c))  
	isADMG(a1)    # Not an AG. (a2) p.969    
	a2 <- makeMG(dg=DAG(b ~ a, d~c), bg=UG(~a*c+c*b+b*d))           # Fig. 3 (b1) p.969  
	isADMG(a2)



cleanEx()
nameEx("isAG")
### * isAG

flush(stderr()); flush(stdout())

### Name: isAG
### Title: Ancestral graph
### Aliases: isAG
### Keywords: graphs ancestral graph mixed graph models multivariate

### ** Examples

	## Examples from Richardson and Spirtes (2002)
	a1 <- makeMG(dg=DAG(a~b, b~d, d~c), bg=UG(~a*c))  
	isAG(a1)    # Not an AG. (a2) p.969    
	a2 <- makeMG(dg=DAG(b ~ a, d~c), bg=UG(~a*c+c*b+b*d))           # Fig. 3 (b1) p.969  
	isAG(a2)



cleanEx()
nameEx("isAcyclic")
### * isAcyclic

flush(stderr()); flush(stdout())

### Name: isAcyclic
### Title: Graph queries
### Aliases: isAcyclic
### Keywords: graphs models multivariate

### ** Examples

## A cyclic graph
d <- matrix(0,3,3)
rownames(d) <- colnames(d) <- c("x", "y", "z")
d["x","y"] <- d["y", "z"] <- d["z", "x"] <- 1
## Test if the graph is acyclic
isAcyclic(d)



cleanEx()
nameEx("isGident")
### * isGident

flush(stderr()); flush(stdout())

### Name: isGident
### Title: G-identifiability of an UG
### Aliases: isGident
### Keywords: graphs models multivariate

### ** Examples

## A not G-identifiable UG
G1 <- UG(~ a*b + u*v)
isGident(G1)
## G-identifiable UG
G2 <- UG(~ a + b + u*v)
isGident(G2)
## G-identifiable UG
G3 <- cmpGraph(UG(~a*b*c+x*y*z))
isGident(G3)



cleanEx()
nameEx("makeMG")
### * makeMG

flush(stderr()); flush(stdout())

### Name: makeMG
### Title: Mixed Graphs
### Aliases: makeMG
### Keywords: graphs ancestral graph mixed graph models multivariate

### ** Examples

## Examples from Richardson and Spirtes (2002)
a1 <- makeMG(dg=DAG(a~b, b~d, d~c), bg=UG(~a*c))  
isAG(a1)    # Not an AG. (a2) p.969    
a2 <- makeMG(dg=DAG(b ~ a, d~c), bg=UG(~a*c+c*b+b*d))           # Fig. 3 (b1) p.969  
isAG(a1)
a3 <- makeMG(ug = UG(~ a*c), dg=DAG(b ~ a, d~c), bg=UG(~ b*d)) # Fig. 3 (b2) p.969
a5 <- makeMG(bg=UG(~alpha*beta+gamma*delta), dg=DAG(alpha~gamma,
delta~beta))  # Fig. 6 p. 973
## Another Example
a4 <- makeMG(ug=UG(~y0*y1), dg=DAG(y4~y2, y2~y1), bg=UG(~y2*y3+y3*y4))  
## A mixed graphs with double edges. 
mg <- makeMG(dg = DG(Y ~ X, Z~W, W~Z, Q~X), ug = UG(~X*Q), bg = UG(~ Y*X+X*Q+Q*W + Y*Z) )
## Chronic pain data: a regression graph
chronic.pain <- makeMG(dg = DAG(Y ~ Za, Za ~ Zb + A, Xa ~ Xb, Xb ~ U+V, U ~ A + V, Zb ~ B, A ~ B), bg = UG(~Za*Xa + Zb*Xb))



cleanEx()
nameEx("marg.param")
### * marg.param

flush(stderr()); flush(stdout())

### Name: marg.param
### Title: Link function of marginal log-linear parameterization
### Aliases: marg.param
### Keywords: logistic models ordinal models

### ** Examples
    
marg.param(c(3,3), c("l", "g"))



cleanEx()
nameEx("marks")
### * marks

flush(stderr()); flush(stdout())

### Name: marks
### Title: Mathematics marks
### Aliases: marks
### Keywords: datasets

### ** Examples

data(marks)
pairs(marks)



cleanEx()
nameEx("mat.mlogit")
### * mat.mlogit

flush(stderr()); flush(stdout())

### Name: mat.mlogit
### Title: Multivariate logistic parametrization
### Aliases: mat.mlogit
### Keywords: logistic model

### ** Examples
 
mat.mlogit(2)



cleanEx()
nameEx("msep")
### * msep

flush(stderr()); flush(stdout())

### Name: msep
### Title: The m-separation criterion
### Aliases: msep
### Keywords: graphs d-separation m-separation mixed graph

### ** Examples

H <-matrix(c(0,0,0,0,
	         1,0,0,1,
	         0,1,0,0,
	         0,0,0,0),4,4)
msep(H,1,4, 2)
msep(H,1,4, c())



cleanEx()
nameEx("null")
### * null

flush(stderr()); flush(stdout())

### Name: null
### Title: Null space of a matrix
### Aliases: null
### Keywords: matrix

### ** Examples

 null(c(1,1,1))



cleanEx()
nameEx("parcor")
### * parcor

flush(stderr()); flush(stdout())

### Name: parcor
### Title: Partial correlations
### Aliases: parcor
### Keywords: array graphs models multivariate

### ** Examples

### Partial correlations for the mathematics marks data
data(marks)
S <- var(marks)
parcor(S)



cleanEx()
nameEx("pcor")
### * pcor

flush(stderr()); flush(stdout())

### Name: pcor
### Title: Partial correlation
### Aliases: pcor
### Keywords: models multivariate

### ** Examples

data(marks)
## The correlation between vectors and algebra given analysis and statistics
 pcor(c("vectors", "algebra", "analysis", "statistics"), var(marks))
## The same
pcor(c(2,3,4,5), var(marks))
## The correlation between vectors and algebra given statistics
 pcor(c("vectors", "algebra", "statistics"), var(marks))
## The marginal correlation between analysis and statistics 
pcor(c("analysis","statistics"), var(marks))



cleanEx()
nameEx("pcor.test")
### * pcor.test

flush(stderr()); flush(stdout())

### Name: pcor.test
### Title: Test for zero partial association
### Aliases: pcor.test
### Keywords: htest multivariate

### ** Examples

## Are 2,3 independent given 1?
data(marks)
pcor.test(pcor(c(2,3,1), var(marks)), 1, n=88)



cleanEx()
nameEx("plotGraph")
### * plotGraph

flush(stderr()); flush(stdout())

### Name: plotGraph
### Title: Plot of a mixed graph
### Aliases: plotGraph
### Keywords: graphs adjacency matrix mixed graphs plot

### ** Examples

exvec<-c("b",1,2,"b",1,14,"a",9,8,"l",9,11,
         "a",10,8,"a",11,2,"a",11,9,"a",11,10,
         "a",12,1,"b",12,14,"a",13,10,"a",13,12)
plotGraph(exvec)
############################################
amat<-matrix(c(0,11,0,0,10,0,100,0,0,100,0,1,0,0,1,0),4,4)
plotGraph(amat)     
plotGraph(makeMG(bg = UG(~a*b*c+ c*d), dg = DAG(a ~ x + z, b ~ z )))
plotGraph(makeMG(bg = UG(~a*b*c+ c*d), dg = DAG(a ~ x + z, b ~ z )), dashed = TRUE)    
# A graph with double and triple edges
G <-
structure(c(0, 101, 0, 0, 100, 0, 100, 100, 0, 100, 0, 100, 0, 
111, 100, 0), .Dim = c(4L, 4L), .Dimnames = list(c("X", "Z", 
"Y", "W"), c("X", "Z", "Y", "W")))
plotGraph(G)      
# A regression chain graph with longer labels
 plotGraph(makeMG(bg = UG(~Love*Constraints+ Constraints*Reversal+ Abuse*Distress), 
   dg = DAG(Love ~ Abuse + Distress, Constraints ~ Distress, Reversal ~ Distress, 
   Abuse ~ Fstatus, Distress ~ Fstatus), 
   ug = UG(~Fstatus*Schooling+ Schooling*Age)), 
   dashed = TRUE, noframe = TRUE)    
# A graph with 4 edges between two nodes. 
G4 = matrix(0, 2, 2); G4[1,2] = 111; G4[2,1] = 111
plotGraph(G4)



cleanEx()
nameEx("powerset")
### * powerset

flush(stderr()); flush(stdout())

### Name: powerset
### Title: Power set
### Aliases: powerset
### Keywords: sets

### ** Examples

powerset(c("A", "B", "C"), nonempty = FALSE)  
powerset(1:3, sort = FALSE, nonempty = TRUE)



cleanEx()
nameEx("rcorr")
### * rcorr

flush(stderr()); flush(stdout())

### Name: rcorr
### Title: Random correlation matrix
### Aliases: rcorr
### Keywords: distribution multivariate

### ** Examples

## A random correlation matrix of order 3
rcorr(3)
## A random correlation matrix of order 5
rcorr(5)



cleanEx()
nameEx("rnormDag")
### * rnormDag

flush(stderr()); flush(stdout())

### Name: rnormDag
### Title: Random sample from a decomposable Gaussian model
### Aliases: rnormDag
### Keywords: distribution multivariate

### ** Examples

## Generate a sample of 100 observation from a multivariate normal
## The matrix of the path coefficients 
A <- matrix(
c(1, -2, -3,  0, 0,  0,  0,
  0,  1,  0, -4, 0,  0,  0,
  0,  0,  1,  2, 0,  0,  0,
  0,  0,  0,  1, 1, -5,  0,
  0,  0,  0,  0, 1,  0,  3,
  0,  0,  0,  0, 0,  1, -4,
  0,  0,  0,  0, 0,  0,  1), 7, 7, byrow=TRUE)
D <- rep(1, 7)
X <- rnormDag(100, A, D)

## The true covariance matrix
solve(A) %*% diag(D) %*% t(solve(A))

## Triangular decomposition of the sample covariance matrix
triDec(cov(X))$A



cleanEx()
nameEx("rsphere")
### * rsphere

flush(stderr()); flush(stdout())

### Name: rsphere
### Title: Random vectors on a sphere
### Aliases: rsphere
### Keywords: distribution multivariate

### ** Examples

## 100 points on circle
z <- rsphere(100,2)
plot(z)

## 100 points on a sphere
z <- rsphere(100, 3)
pairs(z)



cleanEx()
nameEx("shipley.test")
### * shipley.test

flush(stderr()); flush(stdout())

### Name: shipley.test
### Title: Test of all independencies implied by a given DAG
### Aliases: shipley.test
### Keywords: graphs models multivariate

### ** Examples

## A decomposable model for the mathematics marks data
data(marks)
dag <- DAG(mechanics ~ vectors+algebra, vectors ~ algebra, statistics ~ algebra+analysis, analysis ~ algebra)
shipley.test(dag, cov(marks), n=88)



cleanEx()
nameEx("stress")
### * stress

flush(stderr()); flush(stdout())

### Name: stress
### Title: Stress
### Aliases: stress
### Keywords: datasets

### ** Examples

data(stress)
G = UG(~ Y*X + X*V + V*U + U*Y)
fitConGraph(G, stress, 100)



cleanEx()
nameEx("surdata")
### * surdata

flush(stderr()); flush(stdout())

### Name: surdata
### Title: A simulated data set
### Aliases: surdata
### Keywords: datasets

### ** Examples

data(surdata)
pairs(surdata)



cleanEx()
nameEx("swp")
### * swp

flush(stderr()); flush(stdout())

### Name: swp
### Title: Sweep operator
### Aliases: swp
### Keywords: array algebra models multivariate

### ** Examples

## A very simple example
V <- matrix(c(10, 1, 1, 2), 2, 2)
swp(V, 2)



cleanEx()
nameEx("topSort")
### * topSort

flush(stderr()); flush(stdout())

### Name: topSort
### Title: Topological sort
### Aliases: topSort topOrder
### Keywords: graphs models multivariate

### ** Examples

## A simple example
dag <- DAG(a ~ b, c ~ a + b, d ~ c + b)
dag
topOrder(dag)
topSort(dag)



cleanEx()
nameEx("transClos")
### * transClos

flush(stderr()); flush(stdout())

### Name: transClos
### Title: Transitive closure of a graph
### Aliases: transClos
### Keywords: graphs models multivariate

### ** Examples

## Closure of a DAG
d <- DAG(y ~ x, x ~ z)
transClos(d)

## Closure of an UG
g <- UG(~ x*y*z+z*u+u*v)
transClos(g)



cleanEx()
nameEx("triDec")
### * triDec

flush(stderr()); flush(stdout())

### Name: triDec
### Title: Triangular decomposition of a covariance matrix
### Aliases: triDec
### Keywords: array algebra models multivariate

### ** Examples

## Triangular decomposition of a covariance matrix
B <- matrix(c(1,  -2, 0, 1,
              0,   1, 0, 1,
              0,   0, 1, 0,
              0,   0, 0, 1), 4, 4, byrow=TRUE)
B
D <- diag(c(3, 1, 2, 1))
S <- B %*% D %*% t(B)
triDec(S)
solve(B)



cleanEx()
nameEx("unmakeMG")
### * unmakeMG

flush(stderr()); flush(stdout())

### Name: unmakeMG
### Title: Loopless mixed graphs components
### Aliases: unmakeMG
### Keywords: graphs ancestral graph mixed graph models multivariate

### ** Examples

ag <- makeMG(ug=UG(~y0*y1), dg=DAG(y4~y2, y2~y1), bg=UG(~y2*y3+y3*y4))  
isAG(ag)
unmakeMG(ag)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
