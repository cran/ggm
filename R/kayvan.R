#### Functions by Kayvan Sadeghi 2011-2012
## May 2012 Changed the return values of some functions to TRUE FALSE
require(graph)
require(igraph)

######################################################################
######################################################################
rem<-function(a,r){ # this is setdiff (a, r)
	k<-0
	b<-a
	for (i in a){
		k<-k+1
		for(j in r){
			if(i==j){
				b<-b[-k]
				k<-k-1
				break}

			}
		}
	return(b)
}
#############################################################################
#############################################################################
#SPl<-function(a,alpha){
#	a<-sort(a)
#	alpha<-sort(alpha)
#	r<-c()
#	if (length(alpha)>0){
#		for(i in 1:length(a)){
#			for(j in 1:length(alpha)){
#				if(a[i]==alpha[j]){
#					r<-c(r,i)
#					break}}}}
#	return(r)
#}
###############################################################################
# Finds indices of b in sorted a

'SPl' = function(a, b) (1:length(a))[is.element(sort(a), b)]
##############################################################################
RR<-function(a){   ## This is unique(a)
	a<-sort(a)
	r<-a[1]
	i<-1
	while(i<length(a)){
		if(a[i]==a[i+1] ){
			i<-i+1}
		else{
			r<-c(r,a[i+1])
			i<-i+1
			}}
	return(r)
}
#################################################################################
#################################################################################


#' Graph to adjacency matrix
#'
#' \code{grMAT} generates the associated adjacency matrix to a given graph.
#'
#'
#' @param agr A graph that can be a \code{graphNEL} or an \code{\link{igraph}}
#' object or a vector of length \eqn{3e}, where \eqn{e} is the number of edges
#' of the graph, that is a sequence of triples (type, node1label, node2label).
#' The type of edge can be \code{"a"} (arrows from node1 to node2), \code{"b"}
#' (arcs), and \code{"l"} (lines).
#' @return A matrix that consists 4 different integers as an \eqn{ij}-element:
#' 0 for a missing edge between \eqn{i} and \eqn{j}, 1 for an arrow from
#' \eqn{i} to \eqn{j}, 10 for a full line between \eqn{i} and \eqn{j}, and 100
#' for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers are added
#' to be associated with multiple edges of different types. The matrix is
#' symmetric w.r.t full lines and bi-directed arrows.
#' @author Kayvan Sadeghi
#' @keywords graphs adjacency matrix mixed graph vector
#' @examples
#'
#' ## Generating the adjacency matrix from a vector
#' exvec <-c ('b',1,2,'b',1,14,'a',9,8,'l',9,11,'a',10,8,
#'            'a',11,2,'a',11,10,'a',12,1,'b',12,14,'a',13,10,'a',13,12)
#' grMAT(exvec)
#'
#'
`grMAT` <- function(agr)
{
	 if (class(agr)[1] == "graphNEL") {
        	agr<-igraph.from.graphNEL(agr)
       }
	 if (class(agr)[1]== "igraph"){
		return(get.adjacency(agr, sparse = FALSE))
	 }
	 if(class(agr)[1] == "character"){
		if (length(agr)%%3!=0){
			stop("'The character object' is not in a valid form")}
		seqt<-seq(1,length(agr),3)
		b<-agr[seqt]
		agrn<- agr[-seqt]
		bn<-c()
		for(i in 1:length(b)){
			if(b[i]!="a" && b[i]!="l" & b[i]!="b"){
				stop("'The numeric object' is not in a valid form")}
			if(b[i]=="l"){
				bn[i]<-10}
			if(b[i]=="a"){
				bn[i]<-1}
			if(b[i]=="b"){
				bn[i]<-100}
			}
		Ragr<-RR(agrn)
		ma<-length(Ragr)
		mat<-matrix(rep(0,(ma)^2),ma,ma)
		for(i in seq(1,length(agrn),2)){
			if((bn[(i+1)/2]==1 & mat[SPl(Ragr,agrn[i]),SPl(Ragr,agrn[i+1])]%%10!=1)|(bn[(i+1)/2]==10 & mat[SPl(Ragr,agrn[i]),SPl(Ragr,agrn[i+1])]%%100<10)|(bn[(i+1)/2]==100 & mat[SPl(Ragr,agrn[i]),SPl(Ragr,agrn[i+1])]<100)){
				mat[SPl(Ragr,agrn[i]),SPl(Ragr,agrn[i+1])]<-mat[SPl(Ragr,agrn[i]),SPl(Ragr,agrn[i+1])]+bn[(i+1)/2]
			if(bn[(i+1)/2]==10 | bn[(i+1)/2]==100){
				mat[SPl(Ragr,agrn[i+1]),SPl(Ragr,agrn[i])]<-mat[SPl(Ragr,agrn[i+1]),SPl(Ragr,agrn[i])]+bn[(i+1)/2]}}
		}
		rownames(mat)<-Ragr
		colnames(mat)<-Ragr
	}
	return(mat)
}


#' Ribbonless graph
#'
#' \code{RG} generates and plots ribbonless graphs (a modification of MC graph
#' to use m-separation) after marginalization and conditioning.
#'
#'
#' @param amat An adjacency matrix, or a graph that can be a \code{graphNEL} or
#' an \code{\link{igraph}} object or a vector of length \eqn{3e}, where \eqn{e}
#' is the number of edges of the graph, that is a sequence of triples (type,
#' node1label, node2label). The type of edge can be \code{"a"} (arrows from
#' node1 to node2), \code{"b"} (arcs), and \code{"l"} (lines).
#' @param M A subset of the node set of \code{a} that is going to be
#' marginalized over
#' @param C Another disjoint subset of the node set of \code{a} that is going
#' to be conditioned on.
#' @param showmat A logical value. \code{TRUE} (by default) to print the
#' generated matrix.
#' @param plot A logical value, \code{FALSE} (by default). \code{TRUE} to plot
#' the generated graph.
#' @param plotfun Function to plot the graph when \code{plot == TRUE}. Can be
#' \code{plotGraph} (the default) or \code{drawGraph}.
#' @param \dots Further arguments passed to \code{plotfun}.
#' @return A matrix that consists 4 different integers as an \eqn{ij}-element:
#' 0 for a missing edge between \eqn{i} and \eqn{j}, 1 for an arrow from
#' \eqn{i} to \eqn{j}, 10 for a full line between \eqn{i} and \eqn{j}, and 100
#' for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers are added
#' to be associated with multiple edges of different types. The matrix is
#' symmetric w.r.t full lines and bi-directed arrows.
#' @author Kayvan Sadeghi
#' @seealso \code{\link{AG}},, \code{\link{MRG}}, \code{\link{SG}}
#' @references Koster, J.T.A. (2002). Marginalizing and conditioning in
#' graphical models.  \emph{Bernoulli}, 8(6), 817-840.
#'
#' Sadeghi, K. (2011). Stable classes of graphs containing directed acyclic
#' graphs.  \emph{Submitted}.
#' @keywords graphs directed acyclic graph marginalisation and conditioning MC
#' graph ribbonless graph
#' @examples
#'
#' 	ex <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ##The adjacency matrix of a DAG
#' 	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
#' 	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
#' 	               0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
#' 	               0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,
#' 	               1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,
#' 	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
#' 	               0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0),16,16, byrow = TRUE)
#'
#' M<-c(3,5,6,15,16)
#' C<-c(4,7)
#' RG(ex,M,C,plot=TRUE)
#'
`RG` <- function (amat,M=c(),C=c(),showmat=TRUE,plot=FALSE, plotfun = plotGraph, ...)
{
	if(class(amat)[1] == "igraph" || class(amat)[1] == "graphNEL" || class(amat)[1] == "character") {
		amat<-grMAT(amat)}
	if(is(amat,"matrix")){
		if(nrow(amat)==ncol(amat)){
			if(length(rownames(amat))!=ncol(amat)){
	 			rownames(amat)<-1:ncol(amat)
	 			colnames(amat)<-1:ncol(amat)}
					}
		else {
      	  stop("'object' is not in a valid adjacency matrix form")}}
	if(!is(amat,"matrix")) {
	stop("'object' is not in a valid form")}


	S<-C
	St<-c()
	while(identical(S,St)==FALSE)
	{
		St<-S
		for(j in S){
			for(i in rownames(amat)){
				if(amat[i,j]%%10 == 1){
					S<-c(i,S[S!=i])}
				}
			}
	}


	amatr<-amat
	amatt<-2*amat
	while(identical(amatr,amatt)==FALSE)
	{
		amatt<-amatr

############################################################1
			amat21 <- amatr
 	   for (kk in M) {
		  for (kkk in 1:ncol(amatr)) {
			amat31<-amatr
		 	if (amatr[kkk,kk]%%10==1|amatr[kk,kkk]%%10==1) {
				amat31[kk,kkk]<-amatr[kkk,kk]
				amat31[kkk,kk]<-amatr[kk,kkk]}

       		 idx <- which(amat31[, kk]%%10 == 1)
        		lenidx <- length(idx)
        		if ((lenidx > 1&amat31[kkk,kk]%%10==1)) {
           	 	for (ii in 1:(lenidx )) {
               		 #for (jj in (ii + 1):lenidx) {
                 	 	if(amatr[idx[ii], kkk]%%10==0&idx[ii]!=kkk){
			 	amat21[idx[ii], kkk] <- TRUE
                		}}
            		}
       		 }
			}

		amatr<-amat21

################################################################2

		amat22 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in S) {
        		idx <- which(amatr[, kk]%%10 == 1)
			idy <- which(amatr[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%10==0 & idx[ii]!=idy[jj]){
                  amat22[idx[ii], idy[jj]] <- 1
      				          }}
	            		}
        			}
 		   }
		amatr<-amat22+amatr

#################################################################3

		amat33<- t(amatr)
		amat23 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in M) {
        		idx <- which(amat33[, kk]%%10 == 1)
			idy <- which(amat33[, kk]%%100> 9)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idy[jj], idx[ii]]%%10==0 & idx[ii]!=idy[jj]){
                  amat23[idy[jj], idx[ii]] <- 1
      				          }}
	            		}
        			}
 		   }
		amatr<-amat23+amatr

##################################################################4
			amat24 <- matrix(rep(0,length(amat)),dim(amat))
    		 	for (kk in M) {
        		idx <- which(amatr[, kk]%%100>9)
			idy <- which(amatr[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%10==0 & idx[ii]!=idy[jj]){
                  amat24[idx[ii], idy[jj]] <- 1
      				          }}
	            		}
        			}
 		   }
		amatr<-amat24+amatr

####################################################################5



		 amat35<- t(amatr)
   		 amat25 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in M) {
     			   idx <- which(amat35[, kk]%%10 == 1)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]<100){
				   amat25[idx[ii], idx[jj]] <- 100
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat25+t(amat25)+amatr

######################################################################6
		amat36<- t(amatr)
		amat26 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in M) {
        		idx <- which(amat36[, kk]%%10 == 1)
			idy <- which(amat36[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idy[jj], idx[ii]]<100 & idx[ii]!=idy[jj]){
                  amat26[idy[jj], idx[ii]] <- 100
      				          }}
	            		}
        			}
 		   }
		 amatr<-amat26+t(amat26)+amatr

#################################################################7

   		 amat27 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in S) {
     			   idx <- which(amatr[, kk]>99)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]<100){
				   amat27[idx[ii], idx[jj]] <- 100
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat27+t(amat27)+amatr

################################################################8

             amat28 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in S) {
        		idx <- which(amatr[, kk]%%10 == 1)
        		lenidx <- length(idx)
        		if ((lenidx > 1)) {
            	for (ii in 1:(lenidx - 1)) {
                	for (jj in (ii + 1):lenidx) {
			if(amatr[idx[ii], idx[jj]]%%100<10){
                  amat28[idx[ii], idx[jj]] <- 10
      				          }}
	            		}
        			}
 		   }
		amatr<-amat28+t(amat28)+amatr

#################################################################9

			amat29 <- matrix(rep(0,length(amat)),dim(amat))
    		 	for (kk in M) {
        		idx <- which(amatr[, kk]%%10 == 1)
			idy <- which(amatr[, kk]%%100> 9)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%100<10 & idx[ii]!=idy[jj]){
                  amat29[idx[ii], idy[jj]] <- 10
      				          }}
	            		}
        			}
 		   }
		amatr<-amat29+t(amat29)+amatr

##################################################################10

   		 amat20 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in M) {
     			   idx <- which(amatr[, kk]%%100>9)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]%%100<10){
				   amat20[idx[ii], idx[jj]] <- 10
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat20+t(amat20)+amatr

	}
	Mn<-c()
	Cn<-c()
	for(i in M){
		Mn<-c(Mn,which(rownames(amat)==i))}
	for(i in C){
		Cn<-c(Cn,which(rownames(amat)==i))}
	if(length(Mn)>0&length(Cn)>0){
		fr<-amatr[-c(Mn,Cn),-c(Mn,Cn)]}
	if(length(Mn)>0&length(Cn)==0){
		fr<-amatr[-c(Mn),-c(Mn)]}
	if(length(Mn)==0&length(Cn)>0){
		fr<-amatr[-c(Cn),-c(Cn)]}
	if(length(Mn)==0&length(Cn)==0){
		fr<-amatr}
	if(plot==TRUE){
		plotfun(fr,...)}
	if(showmat==FALSE){
		invisible(fr)}
	else{return(fr)}
}
##############################################################################
##############################################################################


#' summary graph
#'
#' \code{SG} generates and plots summary graphs after marginalization and
#' conditioning.
#'
#'
#' @param amat An adjacency matrix, or a graph that can be a \code{graphNEL} or
#' an \code{\link{igraph}} object or a vector of length \eqn{3e}, where \eqn{e}
#' is the number of edges of the graph, that is a sequence of triples (type,
#' node1label, node2label). The type of edge can be \code{"a"} (arrows from
#' node1 to node2), \code{"b"} (arcs), and \code{"l"} (lines).
#' @param M A subset of the node set of \code{a} that is going to be
#' marginalised over
#' @param C Another disjoint subset of the node set of \code{a} that is going
#' to be conditioned on.
#' @param showmat A logical value. \code{TRUE} (by default) to print the
#' generated matrix.
#' @param plot A logical value, \code{FALSE} (by default). \code{TRUE} to plot
#' the generated graph.
#' @param plotfun Function to plot the graph when \code{plot == TRUE}. Can be
#' \code{plotGraph} (the default) or \code{drawGraph}.
#' @param \dots Further arguments passed to \code{plotfun}.
#' @return A matrix that consists 4 different integers as an \eqn{ij}-element:
#' 0 for a missing edge between \eqn{i} and \eqn{j}, 1 for an arrow from
#' \eqn{i} to \eqn{j}, 10 for a full line between \eqn{i} and \eqn{j}, and 100
#' for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers are added
#' to be associated with multiple edges of different types. The matrix is
#' symmetric w.r.t full lines and bi-directed arrows.
#' @author Kayvan Sadeghi
#' @seealso \code{\link{AG}}, \code{\link{MSG}}, \code{\link{RG}}
#' @references Sadeghi, K. (2011). Stable classes of graphs containing directed
#' acyclic graphs.  \emph{Submitted}.
#'
#' Wermuth, N. (2011). Probability distributions with summary graph structure.
#' \emph{Bernoulli}, 17(3),845-879.
#' @keywords graphs directed acyclic graph marginalization and conditioning
#' summary graph
#' @examples
#'
#' 	ex <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ##The adjacency matrix of a DAG
#' 	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
#' 	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
#' 	               0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
#' 	               0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,
#' 	               1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,
#' 	               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#' 	               1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
#' 	               0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0),16,16, byrow = TRUE)
#' M <- c(3,5,6,15,16)
#' C <- c(4,7)
#' SG(ex, M, C, plot = TRUE)
#' SG(ex, M, C, plot = TRUE, plotfun = drawGraph, adjust = FALSE)
#'
`SG`<-function (amat,M=c(),C=c(),showmat=TRUE, plot=FALSE, plotfun = plotGraph, ...)
{
	if(class(amat)[1] == "igraph" || class(amat)[1] == "graphNEL" || class(amat)[1] == "character") {
		amat<-grMAT(amat)}
	if(is(amat,"matrix")){
		if(nrow(amat)==ncol(amat)){
			if(length(rownames(amat))!=ncol(amat)){
	 			rownames(amat)<-1:ncol(amat)
	 			colnames(amat)<-1:ncol(amat)}
					}
		else {
      	  stop("'object' is not in a valid adjacency matrix form")}}
	if(!is(amat,"matrix")) {
	stop("'object' is not in a valid form")}

	S<-C
	St<-c()
	while(identical(S,St)==FALSE)
	{
		St<-S
		for(j in S){
			for(i in rownames(amat)){
				if(amat[i,j]%%10 == 1){
					S<-c(i,S[S!=i])}
				}
			}
	}

	amatr<-amat
	amatt<-2*amat
	while(identical(amatr,amatt)==FALSE)
	{
		amatt<-amatr

############################################################1
			amat21 <- amatr
 	   for (kk in M) {
		  for (kkk in 1:ncol(amatr)) {
			amat31<-amatr
		 	if (amatr[kkk,kk]%%10==1|amatr[kk,kkk]%%10==1) {
				amat31[kk,kkk]<-amatr[kkk,kk]
				amat31[kkk,kk]<-amatr[kk,kkk]}

       		 idx <- which(amat31[, kk]%%10 == 1)
        		lenidx <- length(idx)
        		if ((lenidx > 1&amat31[kkk,kk]%%10==1)) {
           	 	for (ii in 1:(lenidx )) {
               		 #for (jj in (ii + 1):lenidx) {
                 	 	if(amatr[idx[ii], kkk]%%10==0&idx[ii]!=kkk){
			 	amat21[idx[ii], kkk] <- TRUE
                		}}
            		}
       		 }
			}

		amatr<-amat21

################################################################2

		amat22 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in S) {
        		idx <- which(amatr[, kk]%%10 == 1)
			idy <- which(amatr[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%10==0 & idx[ii]!=idy[jj]){
                  amat22[idx[ii], idy[jj]] <- 1
      				          }}
	            		}
        			}
 		   }
		amatr<-amat22+amatr

#################################################################3

		amat33<- t(amatr)
		amat23 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in M) {
        		idx <- which(amat33[, kk]%%10 == 1)
			idy <- which(amat33[, kk]%%100> 9)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idy[jj], idx[ii]]%%10==0 & idx[ii]!=idy[jj]){
                  amat23[idy[jj], idx[ii]] <- 1
      				          }}
	            		}
        			}
 		   }
		amatr<-amat23+amatr

##################################################################4
			amat24 <- matrix(rep(0,length(amat)),dim(amat))
    		 	for (kk in M) {
        		idx <- which(amatr[, kk]%%100>9)
			idy <- which(amatr[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%10==0 & idx[ii]!=idy[jj]){
                  amat24[idx[ii], idy[jj]] <- 1
      				          }}
	            		}
        			}
 		   }
		amatr<-amat24+amatr

####################################################################5



		 amat35<- t(amatr)
   		 amat25 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in M) {
     			   idx <- which(amat35[, kk]%%10 == 1)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]<100){
				   amat25[idx[ii], idx[jj]] <- 100
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat25+t(amat25)+amatr

######################################################################6
		amat36<- t(amatr)
		amat26 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in M) {
        		idx <- which(amat36[, kk]%%10 == 1)
			idy <- which(amat36[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idy[jj], idx[ii]]<100 & idx[ii]!=idy[jj]){
                  amat26[idy[jj], idx[ii]] <- 100
      				          }}
	            		}
        			}
 		   }
		 amatr<-amat26+t(amat26)+amatr

#################################################################7

   		 amat27 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in S) {
     			   idx <- which(amatr[, kk]>99)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]<100){
				   amat27[idx[ii], idx[jj]] <- 100
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat27+t(amat27)+amatr

################################################################8

             amat28 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in S) {
        		idx <- which(amatr[, kk]%%10 == 1)
        		lenidx <- length(idx)
        		if ((lenidx > 1)) {
            	for (ii in 1:(lenidx - 1)) {
                	for (jj in (ii + 1):lenidx) {
			if(amatr[idx[ii], idx[jj]]%%100<10){
                  amat28[idx[ii], idx[jj]] <- 10
      				          }}
	            		}
        			}
 		   }
		amatr<-amat28+t(amat28)+amatr

#################################################################9

			amat29 <- matrix(rep(0,length(amat)),dim(amat))
    		 	for (kk in M) {
        		idx <- which(amatr[, kk]%%10 == 1)
			idy <- which(amatr[, kk]%%100> 9)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%100<10 & idx[ii]!=idy[jj]){
                  amat29[idx[ii], idy[jj]] <- 10
      				          }}
	            		}
        			}
 		   }
		amatr<-amat29+t(amat29)+amatr

##################################################################10

   		 amat20 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in M) {
     			   idx <- which(amatr[, kk]%%100>9)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]%%100<10){
				   amat20[idx[ii], idx[jj]] <- 10
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat20+t(amat20)+amatr

	}

		for(i in 1:ncol(amatr)) {
		for(j in 1:ncol(amatr)) {
			if(amatr[i,j]%%100>9){
				amatr[i,j]<-10
				for(k in 1:ncol(amatr)){
					if(amatr[k,j]==100){
						amatr[j,k]<-1
						amatr[k,j]<-0
						}}
				}
			}}

	SS<-S
	SSt<-c()
	while(identical(SS,SSt)==FALSE)
	{
		SSt<-SS
		for(j in SS){
			for(i in rownames(amat)) {
				if(amatr[i,j]%%10 == 1){
					SS<-c(i,SS[SS!=i])}
				}
			}
	}

		for(i in SS){
		for(j in SS) {
			if(amatr[i,j]%%10==1){
				amatr[i,j]<-10
				amatr[j,i]<-10}}}



	Mn<-c()
	Cn<-c()
	for(i in M){
		Mn<-c(Mn,which(rownames(amat)==i))}
	for(i in C){
		Cn<-c(Cn,which(rownames(amat)==i))}
	if(length(Mn)>0&length(Cn)>0){
		fr<-amatr[-c(Mn,Cn),-c(Mn,Cn)]}
	if(length(Mn)>0&length(Cn)==0){
		fr<-amatr[-c(Mn),-c(Mn)]}
	if(length(Mn)==0&length(Cn)>0){
		fr<-amatr[-c(Cn),-c(Cn)]}
	if(length(Mn)==0&length(Cn)==0){
		fr<-amatr}
	if(plot==TRUE){
		plotfun(fr,...)}
	if(showmat==FALSE){
		invisible(fr)}
	else{return(fr)}

}
##############################################################################
##############################################################################


#' Ancestral graph
#'
#' \code{AG} generates and plots ancestral graphs after marginalization and
#' conditioning.
#'
#'
#' @param amat An adjacency matrix, or a graph that can be of class
#' \code{graphNEL-class} or an \code{\link{igraph}} object, or a vector of
#' length \eqn{3e}, where \eqn{e} is the number of edges of the graph, that is
#' a sequence of triples (type, node1label, node2label). The type of edge can
#' be \code{"a"} (arrows from node1 to node2), \code{"b"} (arcs), and
#' \code{"l"} (lines).
#' @param M A subset of the node set of \code{a} that is going to be
#' marginalized over
#' @param C Another disjoint subset of the node set of \code{a} that is going
#' to be conditioned on.
#' @param showmat A logical value. \code{TRUE} (by default) to print the
#' generated matrix.
#' @param plot A logical value, \code{FALSE} (by default). \code{TRUE} to plot
#' the generated graph.
#' @param plotfun Function to plot the graph when \code{plot == TRUE}. Can be
#' \code{plotGraph} (the default) or \code{drawGraph}.
#' @param \dots Further arguments passed to \code{plotfun}.
#' @return A matrix that is the adjacency matrix of the generated graph. It
#' consists of 4 different integers as an \eqn{ij}-element: 0 for a missing
#' edge between \eqn{i} and \eqn{j}, 1 for an arrow from \eqn{i} to \eqn{j}, 10
#' for a full line between \eqn{i} and \eqn{j}, and 100 for a bi-directed arrow
#' between \eqn{i} and \eqn{j}. These numbers are added to be associated with
#' multiple edges of different types. The matrix is symmetric w.r.t full lines
#' and bi-directed arrows.
#' @author Kayvan Sadeghi
#' @seealso \code{\link{MAG}}, \code{\link{RG}}, \code{\link{SG}}
#' @references Richardson, T.S. and Spirtes, P. (2002).  Ancestral graph Markov
#' models. \emph{Annals of Statistics}, 30(4), 962-1030.
#'
#' Sadeghi, K. (2011).  Stable classes of graphs containing directed acyclic
#' graphs.  \emph{Submitted}.
#' @keywords graphs ancestral graph directed acyclic graph marginalization and
#' conditioning
#' @examples
#'
#' ##The adjacency matrix of a DAG
#' ex<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
#'              0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,
#'              1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0),16,16,byrow=TRUE)
#' M <- c(3,5,6,15,16)
#' C <- c(4,7)
#' AG(ex, M, C, plot = TRUE)
#'
`AG`<-function (amat,M=c(),C=c(),showmat=TRUE,plot=FALSE, plotfun = plotGraph, ...)
{
	if(class(amat)[1] == "igraph" || class(amat)[1] == "graphNEL" || class(amat)[1] == "character") {
		amat<-grMAT(amat)}
	if(is(amat,"matrix")){
		if(nrow(amat)==ncol(amat)){
			if(length(rownames(amat))!=ncol(amat)){
	 			rownames(amat)<-1:ncol(amat)
	 			colnames(amat)<-1:ncol(amat)}
					}
		else {
      	  stop("'object' is not in a valid adjacency matrix form")}}
	if(!is(amat,"matrix")) {
	stop("'object' is not in a valid form")}

	S<-C
	St<-c()
	while(identical(S,St)==FALSE)
	{
		St<-S
		for(j in S){
			for(i in rownames(amat)){
				if(amat[i,j]%%10 == 1){
					S<-c(i,S[S!=i])}
				}
			}
	}

	amatr<-amat
	amatt<-2*amat
	while(identical(amatr,amatt)==FALSE)
	{
		amatt<-amatr


############################################################1
			amat21 <- amatr
 	   for (kk in M) {
		  for (kkk in 1:ncol(amatr)) {
			amat31<-amatr
		 	if (amatr[kkk,kk]%%10==1|amatr[kk,kkk]%%10==1) {
				amat31[kk,kkk]<-amatr[kkk,kk]
				amat31[kkk,kk]<-amatr[kk,kkk]}

       		 idx <- which(amat31[, kk]%%10 == 1)
        		lenidx <- length(idx)
        		if ((lenidx > 1&amat31[kkk,kk]%%10==1)) {
           	 	for (ii in 1:(lenidx )) {
               		 #for (jj in (ii + 1):lenidx) {
                 	 	if(amatr[idx[ii], kkk]%%10==0&idx[ii]!=kkk){
			 	amat21[idx[ii], kkk] <- TRUE
                		}}
            		}
       		 }
			}

		amatr<-amat21

################################################################2

		amat22 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in S) {
        		idx <- which(amatr[, kk]%%10 == 1)
			idy <- which(amatr[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%10==0 & idx[ii]!=idy[jj]){
                  amat22[idx[ii], idy[jj]] <- 1
      				          }}
	            		}
        			}
 		   }
		amatr<-amat22+amatr

#################################################################3

		amat33<- t(amatr)
		amat23 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in M) {
        		idx <- which(amat33[, kk]%%10 == 1)
			idy <- which(amat33[, kk]%%100> 9)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idy[jj], idx[ii]]%%10==0 & idx[ii]!=idy[jj]){
                  amat23[idy[jj], idx[ii]] <- 1
      				          }}
	            		}
        			}
 		   }
		amatr<-amat23+amatr

##################################################################4
			amat24 <- matrix(rep(0,length(amat)),dim(amat))
    		 	for (kk in M) {
        		idx <- which(amatr[, kk]%%100>9)
			idy <- which(amatr[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%10==0 & idx[ii]!=idy[jj]){
                  amat24[idx[ii], idy[jj]] <- 1
      				          }}
	            		}
        			}
 		   }
		amatr<-amat24+amatr

####################################################################5



		 amat35<- t(amatr)
   		 amat25 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in M) {
     			   idx <- which(amat35[, kk]%%10 == 1)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]<100){
				   amat25[idx[ii], idx[jj]] <- 100
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat25+t(amat25)+amatr

######################################################################6
		amat36<- t(amatr)
		amat26 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in M) {
        		idx <- which(amat36[, kk]%%10 == 1)
			idy <- which(amat36[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idy[jj], idx[ii]]<100 & idx[ii]!=idy[jj]){
                  amat26[idy[jj], idx[ii]] <- 100
      				          }}
	            		}
        			}
 		   }
		 amatr<-amat26+t(amat26)+amatr

#################################################################7

   		 amat27 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in S) {
     			   idx <- which(amatr[, kk]>99)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]<100){
				   amat27[idx[ii], idx[jj]] <- 100
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat27+t(amat27)+amatr

################################################################8

             amat28 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in S) {
        		idx <- which(amatr[, kk]%%10 == 1)
        		lenidx <- length(idx)
        		if ((lenidx > 1)) {
            	for (ii in 1:(lenidx - 1)) {
                	for (jj in (ii + 1):lenidx) {
			if(amatr[idx[ii], idx[jj]]%%100<10){
                  amat28[idx[ii], idx[jj]] <- 10
      				          }}
	            		}
        			}
 		   }
		amatr<-amat28+t(amat28)+amatr

#################################################################9

			amat29 <- matrix(rep(0,length(amat)),dim(amat))
    		 	for (kk in M) {
        		idx <- which(amatr[, kk]%%10 == 1)
			idy <- which(amatr[, kk]%%100> 9)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%100<10 & idx[ii]!=idy[jj]){
                  amat29[idx[ii], idy[jj]] <- 10
      				          }}
	            		}
        			}
 		   }
		amatr<-amat29+t(amat29)+amatr

##################################################################10

   		 amat20 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in M) {
     			   idx <- which(amatr[, kk]%%100>9)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]%%100<10){
				   amat20[idx[ii], idx[jj]] <- 10
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat20+t(amat20)+amatr

	}

		for(i in 1:ncol(amatr)) {
		for(j in 1:ncol(amatr)) {
			if(amatr[i,j]%%100>9){
				amatr[i,j]<-10
				for(k in 1:ncol(amatr)){
					if(amatr[k,j]==100){
						amatr[j,k]<-1
						amatr[k,j]<-0
						}}
				}
			}}


	SS<-S
	SSt<-c()
	while(identical(SS,SSt)==FALSE)
	{
		SSt<-SS
		for(j in SS){
			for(i in rownames(amat)) {
				if(amatr[i,j]%%10 == 1){
					SS<-c(i,SS[SS!=i])}
				}
			}
	}

		for(i in SS){
		for(j in SS) {
			if(amatr[i,j]%%10==1){
				amatr[i,j]<-10
				amatr[j,i]<-10}}}



		amatn <- amatr
 	   for (kk in 1:ncol(amatr)) {
		  for (kkk in 1:ncol(amatr)) {
			amat3n<-amatr
		 	if (amatr[kkk,kk]%%10==1|amatr[kk,kkk]%%10==1) {
				amat3n[kk,kkk]<-amatr[kkk,kk]
				amat3n[kkk,kk]<-amatr[kk,kkk]}

       		 idx <- which(amat3n[, kk]%%10 == 1)
        		lenidx <- length(idx)
        		if ((lenidx > 1&amat3n[kkk,kk]%%10==1)) {
           	 	for (ii in 1:(lenidx )) {
               		 #for (jj in (ii + 1):lenidx) {
                 	 	if(amatr[idx[ii], kkk]%%10==0&idx[ii]!=kkk){
			 	amatn[idx[ii], kkk] <- TRUE
                		}}
            		}
       		 }
			}


	amatt<-2*amat
	while(identical(amatr,amatt)==FALSE)
	{
		amatt<-amatr


	 amat27n <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in 1:ncol(amatr)) {
     			   idx <- which(amatn[, kk]>99)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]<100 && (amatn[kk,idx[ii]]%%10==1 || amatn[kk,idx[jj]]%%10==1)){
				   amat27n[idx[ii], idx[jj]] <- 100
               					 }}
           					 }
      			  }
   		 }
		amatn<-amat27n+t(amat27n)+amatn
		amatr<-amat27n+t(amat27n)+amatr

		amat22n <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in 1:ncol(amatr)) {
        		idx <- which(amatn[, kk]%%10 == 1)
			idy <- which(amatn[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%10==0 & idx[ii]!=idy[jj]& amatn[kk,idy[jj]]%%10==1){

                  amat22n[idx[ii], idy[jj]] <- 1
      				          }}
	            		}
        			}
 		   }
		amatn<-amat22n+amatn
		amatr<-amat22n+amatr
	}

	for(i in 1:ncol(amatr)) {
		for(j in 1:ncol(amatr)) {
			if(amatr[i,j]==101){
				amatr[i,j]<-1
				amatr[j,i]<-0}}}

	Mn<-c()
	Cn<-c()
	for(i in M){
		Mn<-c(Mn,which(rownames(amat)==i))}
	for(i in C){
		Cn<-c(Cn,which(rownames(amat)==i))}
	if(length(Mn)>0&length(Cn)>0){
		fr<-amatr[-c(Mn,Cn),-c(Mn,Cn)]}
	if(length(Mn)>0&length(Cn)==0){
		fr<-amatr[-c(Mn),-c(Mn)]}
	if(length(Mn)==0&length(Cn)>0){
		fr<-amatr[-c(Cn),-c(Cn)]}
	if(length(Mn)==0&length(Cn)==0){
		fr<-amatr}
	if(plot==TRUE){
		plotfun(fr, ...)}
	if(showmat==FALSE){
		invisible(fr)}
	else{return(fr)}
}

##############################################################################
##############################################################################


#' Maximisation for graphs
#'
#' \code{Max} generates a maximal graph that induces the same independence
#' model from a non-maximal graph.
#'
#' \code{Max} looks for non-adjacent pais of nodes that are connected by
#' primitive inducing paths, and connect such pairs by an appropriate edge.
#'
#' @param amat An adjacency matrix, or a graph that can be a \code{graphNEL} or
#' an \code{\link{igraph}} object or a vector of length \eqn{3e}, where \eqn{e}
#' is the number of edges of the graph, that is a sequence of triples (type,
#' node1label, node2label). The type of edge can be \code{"a"} (arrows from
#' node1 to node2), \code{"b"} (arcs), and \code{"l"} (lines).
#' @return A matrix that consists 4 different integers as an \eqn{ij}-element:
#' 0 for a missing edge between \eqn{i} and \eqn{j}, 1 for an arrow from
#' \eqn{i} to \eqn{j}, 10 for a full line between \eqn{i} and \eqn{j}, and 100
#' for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers are added
#' to be associated with multiple edges of different types. The matrix is
#' symmetric w.r.t full lines and bi-directed arrows.
#' @author Kayvan Sadeghi
#' @seealso \code{\link{MAG}}, \code{\link{MRG}}, \code{\link{msep}},
#' \code{\link{MSG}}
#' @references Richardson, T.S. and Spirtes, P. (2002). Ancestral graph Markov
#' models. \emph{Annals of Statistics}, 30(4), 962-1030.
#'
#' Sadeghi, K. and Lauritzen, S.L. (2011). Markov properties for loopless mixed
#' graphs. \emph{Submitted}. \url{http://arxiv.org/abs/1109.5909}.
#' @keywords graphs loopless mixed graph m-separation maximality
#' @examples
#'
#' H <- matrix(c(  0,100,  1,  0,
#' 	          100,  0,100,  0,
#' 	            0,100,  0,100,
#' 	            0,  1,100,  0), 4, 4)
#' Max(H)
#'
Max<-function(amat)
{
	if(class(amat)[1] == "igraph" || class(amat)[1] == "graphNEL" || class(amat)[1] == "character") {
		amat<-grMAT(amat)}
	if(is(amat,"matrix")){
		if(nrow(amat)==ncol(amat)){
			if(length(rownames(amat))!=ncol(amat)){
	 			rownames(amat)<-1:ncol(amat)
	 			colnames(amat)<-1:ncol(amat)}
					}
		else {
      	  stop("'object' is not in a valid adjacency matrix form")}}
	if(!is(amat,"matrix")) {
	stop("'object' is not in a valid form")}

	na<-ncol(amat)
	at<-which(amat+t(amat)+diag(na)==0,arr.ind=TRUE)
	if(dim(at)[1]>0){
	for(i in 1:dim(at)[1]){
		S<-at[i,]
		St<-c()
		while(identical(S,St)==FALSE){
			St<-S
			for(j in S){
				for(k in 1:na){
					if(amat[k,j]%%10 == 1){
						S<-c(k,S[S!=k])}
					}
				}}
	one<-at[i,1]
	two<-at[i,2]
	onearrow<-c()
	onearc<-c()
	twoarrow<-c()
	twoarc<-c()
	Sr<-S
	Sr<-Sr[Sr!=at[i,1]]
	Sr<-Sr[Sr!=at[i,2]]
	if(length(Sr)>0){
		for(j in Sr){
			if(amat[one,j]%%10==1){
				onearrow<-c(onearrow,j)}
 			if(amat[one,j]>99){
				onearc<-c(onearc,j)}}
		for(j in Sr){
			if(amat[two,j]%%10==1){
				for(k in onearc){
					if(j==k || amat[j,k]>99){
						amat[at[i,2],at[i,1]]<-1}}
				twoarrow<-c(twoarrow,j)}
			if(amat[two,j]>99){
				for(k in onearrow){
					if(j==k || amat[j,k]>99){
						amat[at[i,1],at[i,2]]<-1}}
				for(k in onearc){
					if(j==k || amat[j,k]>99){
						amat[at[i,1],at[i,2]]<-100
						amat[at[i,2],at[i,1]]<-100}}
				twoarc<-c(twoarc,j)}}}
		if(length(c(onearc,onearrow,twoarc,twoarrow))>0){

			for(j in c(onearc,onearrow,twoarc,twoarrow)){
				Sr<-Sr[Sr!=j]}
		onearcn<-c()
		twoarcn<-c()
		onearrown<-c()
		twoarrown<-c()
		while(length(Sr)>0){
			for(l in onearc){
				for(j in Sr){
 					if(amat[l,j]>99){
						onearcn<-c(onearcn,j)}}}
			for(l in onearrow){
				for(j in Sr){
 					if(amat[l,j]>99){
						onearrown<-c(onearrown,j)}}}

			for(l in twoarc){
				for(j in Sr){
					if(amat[l,j]>99){
						for(k in onearrow){
							if(j==k || amat[j,k]>99){
								amat[at[i,1],at[i,2]]<-1}}
						for(k in onearc){
							if(j==k || amat[j,k]>99){
								amat[at[i,1],at[i,2]]<-100
								amat[at[i,2],at[i,1]]<-100}}
						twoarcn<-c(twoarcn,j)}}}
			for(l in twoarrow){
				for(j in Sr){
					if(amat[l,j]>99){
						for(k in onearc){
							if(j==k || amat[j,k]>99){
								amat[at[i,1],at[i,2]]<-100
								amat[at[i,2],at[i,1]]<-100}}
						twoarrown<-c(twoarrown,j)}}}
			if(length(c(onearcn,onearrown,twoarcn,twoarrown))==0){
				break}
			for(j in c(onearcn,onearrown,twoarcn,twoarrown)){
				Sr<-Sr[Sr!=j]}
			onearc<-onearcn
			twoarc<-twoarcn
			onearrow<-onearrown
			twoarrow<-twoarrown
			onearcn<-c()
			twoarcn<-c()
			onearrown<-c()
			twoarrown<-c()}}}}
	return(amat)
}
#####################################################################################################
######################################################################################################


#' The m-separation criterion
#'
#' \code{msep} determines whether two set of nodes are m-separated by a third
#' set of nodes.
#'
#'
#' @param a An adjacency matrix, or a graph that can be a \code{graphNEL} or an
#' \code{\link{igraph}} object or a vector of length \eqn{3e}, where \eqn{e} is
#' the number of edges of the graph, that is a sequence of triples (type,
#' node1label, node2label). The type of edge can be \code{"a"} (arrows from
#' node1 to node2), \code{"b"} (arcs), and \code{"l"} (lines).
#' @param alpha A subset of the node set of \code{a}
#' @param beta Another disjoint subset of the node set of \code{a}
#' @param C A third disjoint subset of the node set of \code{a}
#' @return A logical value. \code{TRUE} if \code{alpha} and \code{beta} are
#' m-separated given \code{C}.  \code{FALSE} otherwise.
#' @author Kayvan Sadeghi
#' @seealso \code{\link{dSep}}, \code{\link{MarkEqMag}}
#' @references Richardson, T.S. and Spirtes, P. (2002) Ancestral graph Markov
#' models. \emph{Annals of Statistics}, 30(4), 962-1030.
#'
#' Sadeghi, K. and Lauritzen, S.L. (2011). Markov properties for loopless mixed
#' graphs. \emph{Submitted}, 2011. URL \url{http://arxiv.org/abs/1109.5909}.
#' @keywords graphs d-separation m-separation mixed graph
#' @examples
#'
#' H <-matrix(c(0,0,0,0,
#' 	         1,0,0,1,
#' 	         0,1,0,0,
#' 	         0,0,0,0),4,4)
#' msep(H,1,4, 2)
#' msep(H,1,4, c())
#'
msep<-function(a,alpha,beta,C=c()){
	if(class(a)[1] == "igraph" || class(a)[1] == "graphNEL" || class(a)[1] == "character") {
		a<-grMAT(a)}
	if(is(a,"matrix")){
		if(nrow(a)==ncol(a)){
			if(length(rownames(a))!=ncol(a)){
	 			rownames(a)<-1:ncol(a)
	 			colnames(a)<-1:ncol(a)}
					}
		else {
      	  stop("'object' is not in a valid adjacency matrix form")}}
	if(!is(a,"matrix")) {
	stop("'object' is not in a valid form")}

		M<-rem(rownames(a),c(alpha,beta,C))
		ar<-Max(RG(a,M,C))

	#aralpha<-as.matrix(ar[SPl(c(alpha,beta),alpha),SPl(c(alpha,beta),beta)])
	#arbeta<-as.matrix(ar[SPl(c(alpha,beta),beta),SPl(c(alpha,beta),alpha)])
	#for(i in 1:length(alpha)){
	#	for(j in 1:length(beta)){
	#		if(aralpha[j,i]!=0 || arbeta[j,i]!=0){
	#			return("NOT separated")
	#			break
	#			break}}}
	if(max(ar[as.character(beta),as.character(alpha)]+ar[as.character(alpha),as.character(beta)]!=0)){
		return(FALSE)}
			return(TRUE)
}
############################################################################
############################################################################


#' Maximal ribbonless graph
#'
#' \code{MRG} generates and plots maximal ribbonless graphs (a modification of
#' MC graph to use m-separation) after marginalisation and conditioning.
#'
#' This function uses the functions \code{\link{RG}} and \code{\link{Max}}.
#'
#' @param amat An adjacency matrix, or a graph that can be a \code{graphNEL} or
#' an \code{\link{igraph}} object or a vector of length \eqn{3e}, where \eqn{e}
#' is the number of edges of the graph, that is a sequence of triples (type,
#' node1label, node2label). The type of edge can be \code{"a"} (arrows from
#' node1 to node2), \code{"b"} (arcs), and \code{"l"} (lines).
#' @param M A subset of the node set of \code{a} that is going to be
#' marginalized over
#' @param C Another disjoint subset of the node set of \code{a} that is going
#' to be conditioned on.
#' @param showmat A logical value. \code{TRUE} (by default) to print the
#' generated matrix.
#' @param plot A logical value, \code{FALSE} (by default). \code{TRUE} to plot
#' the generated graph.
#' @param plotfun Function to plot the graph when \code{plot == TRUE}. Can be
#' \code{plotGraph} (the default) or \code{drawGraph}.
#' @param \dots Further arguments passed to \code{plotfun}.
#' @return A matrix that consists 4 different integers as an \eqn{ij}-element:
#' 0 for a missing edge between \eqn{i} and \eqn{j}, 1 for an arrow from
#' \eqn{i} to \eqn{j}, 10 for a full line between \eqn{i} and \eqn{j}, and 100
#' for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers are added
#' to be associated with multiple edges of different types. The matrix is
#' symmetric w.r.t full lines and bi-directed arrows.
#' @author Kayvan Sadeghi
#' @seealso \code{\link{MAG}}, \code{\link{Max}}, \code{\link{MSG}},
#' \code{\link{RG}}
#' @references Koster, J.T.A. (2002). Marginalizing and conditioning in
#' graphical models.  \emph{Bernoulli}, 8(6), 817-840.
#'
#' Richardson, T.S. and Spirtes, P. (2002). Ancestral graph Markov models.
#' \emph{Annals of Statistics}, 30(4), 962-1030.
#'
#' Sadeghi, K. (2011). Stable classes of graphs containing directed acyclic
#' graphs.  \emph{Submitted}.
#'
#' Sadeghi, K. and Lauritzen, S.L. (2011). Markov properties for loopless mixed
#' graphs. \emph{Submitted}. URL \url{http://arxiv.org/abs/1109.5909}.
#' @keywords graphs directed acyclic graph marginalisation and conditioning
#' maximality of graphs MC graph ribbonless graph
#' @examples
#'
#' ex <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ##The adjacency matrix of a DAG
#'                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'                1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'                0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
#'                0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
#'                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'                0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
#'                0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
#'                0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,
#'                1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'                0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,
#'                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'                1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
#'                0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0),16,16, byrow = TRUE)
#' M <- c(3,5,6,15,16)
#' C <- c(4,7)
#' MRG(ex, M, C, plot = TRUE)
#' ###################################################
#' H <- matrix(c( 0, 100,   1,   0,
#'   	         100,   0, 100,   0,
#'  	             0, 100,   0, 100,
#' 	             0,   1, 100,   0), 4,4)
#' Max(H)
#'
`MRG` <- function (amat,M=c(),C=c(),showmat=TRUE,plot=FALSE, plotfun = plotGraph, ...)
{
	return(Max(RG(amat,M,C,showmat,plot, plotfun = plotGraph, ...)))
}
##########################################################################
##########################################################################


#' Maximal summary graph
#'
#' \code{MAG} generates and plots maximal summary graphs after marginalization
#' and conditioning.
#'
#' This function uses the functions \code{\link{SG}} and \code{\link{Max}}.
#'
#' @param amat An adjacency matrix of a MAG, or a graph that can be a
#' \code{graphNEL} or an \code{\link{igraph}} object or a vector of length
#' \eqn{3e}, where \eqn{e} is the number of edges of the graph, that is a
#' sequence of triples (type, node1label, node2label). The type of edge can be
#' \code{"a"} (arrows from node1 to node2), \code{"b"} (arcs), and \code{"l"}
#' (lines).
#' @param M A subset of the node set of \code{a} that is going to be
#' marginalized over
#' @param C Another disjoint subset of the node set of \code{a} that is going
#' to be conditioned on.
#' @param showmat A logical value. \code{TRUE} (by default) to print the
#' generated matrix.
#' @param plot A logical value, \code{FALSE} (by default). \code{TRUE} to plot
#' the generated graph.
#' @param plotfun Function to plot the graph when \code{plot == TRUE}. Can be
#' \code{plotGraph} (the default) or \code{drawGraph}.
#' @param \dots Further arguments passed to \code{plotfun}.
#' @return A matrix that consists 4 different integers as an \eqn{ij}-element:
#' 0 for a missing edge between \eqn{i} and \eqn{j}, 1 for an arrow from
#' \eqn{i} to \eqn{j}, 10 for a full line between \eqn{i} and \eqn{j}, and 100
#' for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers are added
#' to be associated with multiple edges of different types. The matrix is
#' symmetric w.r.t full lines and bi-directed arrows.
#' @author Kayvan Sadeghi
#' @seealso \code{\link{MAG}}, \code{\link{Max}}, \code{\link{MRG}},
#' \code{\link{SG}}
#' @references Richardson, T.S. and Spirtes, P. (2002). Ancestral graph Markov
#' models. \emph{Annals of Statistics}, 30(4), 962-1030.
#'
#' Sadeghi, K. (2011). Stable classes of graphs containing directed acyclic
#' graphs.  \emph{Submitted}.
#'
#' Sadeghi, K. and Lauritzen, S.L. (2011). Markov properties for loopless mixed
#' graphs. \emph{Submitted}. URL \url{http://arxiv.org/abs/1109.5909}.
#'
#' Wermuth, N. (2011). Probability distributions with summary graph structure.
#' \emph{Bernoulli}, 17(3), 845-879.
#' @keywords graphs directed acyclic graph marginalisation and conditioning
#' maximality of graphs summary graph
#' @examples
#'
#' ex<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ##The adjacency matrix of a DAG
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
#'              0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,
#'              1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0), 16, 16, byrow=TRUE)
#' M <- c(3,5,6,15,16)
#' C <- c(4,7)
#' MSG(ex,M,C,plot=TRUE)
#' ###################################################
#' H<-matrix(c(0,100,1,0,100,0,100,0,0,100,0,100,0,1,100,0),4,4)
#' Max(H)
#'
`MSG` <- function (amat,M=c(),C=c(),showmat=TRUE,plot=FALSE, plotfun = plotGraph, ...)
{
	return(Max(SG(amat,M,C,showmat,plot, plotfun = plotGraph, ...)))
}
############################################################################
###########################################################################


#' Maximal ancestral graph
#'
#' \code{MAG} generates and plots maximal ancestral graphs after
#' marginalisation and conditioning.
#'
#' This function uses the functions \code{\link{AG}} and \code{\link{Max}}.
#'
#' @param amat An adjacency matrix, or a graph that can be a \code{graphNEL} or
#' an \code{\link{igraph}} object or a vector of length \eqn{3e}, where \eqn{e}
#' is the number of edges of the graph, that is a sequence of triples (type,
#' node1label, node2label). The type of edge can be \code{"a"} (arrows from
#' node1 to node2), \code{"b"} (arcs), and \code{"l"} (lines).
#' @param M A subset of the node set of \code{a} that is going to be
#' marginalized over
#' @param C Another disjoint subset of the node set of \code{a} that is going
#' to be conditioned on.
#' @param showmat A logical value. \code{TRUE} (by default) to print the
#' generated matrix.
#' @param plot A logical value, \code{FALSE} (by default). \code{TRUE} to plot
#' the generated graph.
#' @param plotfun Function to plot the graph when \code{plot == TRUE}. Can be
#' \code{plotGraph} (the default) or \code{drawGraph}.
#' @param \dots Further arguments passed to \code{plotfun}.
#' @return A matrix that consists 4 different integers as an \eqn{ij}-element:
#' 0 for a missing edge between \eqn{i} and \eqn{j}, 1 for an arrow from
#' \eqn{i} to \eqn{j}, 10 for a full line between \eqn{i} and \eqn{j}, and 100
#' for a bi-directed arrow between \eqn{i} and \eqn{j}. These numbers are added
#' to be associated with multiple edges of different types. The matrix is
#' symmetric w.r.t full lines and bi-directed arrows.
#' @author Kayvan Sadeghi
#' @seealso \code{\link{AG}}, \code{\link{Max}}, \code{\link{MRG}},
#' \code{\link{MSG}}
#' @references Richardson, T. S. and Spirtes, P. (2002). Ancestral graph Markov
#' models. \emph{Annals of Statistics}, 30(4), 962-1030.
#'
#' Sadeghi, K. (2011). Stable classes of graphs containing directed acyclic
#' graphs.  \emph{Submitted}.
#'
#' Sadeghi, K.  and Lauritzen, S.L. (2011). Markov properties for loopless
#' mixed graphs. \emph{Submitted}. URL \url{http://arxiv.org/abs/1109.5909}.
#' @keywords ancestral graph directed acyclic graph marginalization and
#' conditioning maximality of graphs
#' @examples
#'
#' ex<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ##The adjacency matrix of a DAG
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
#'              0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,
#'              1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'              1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
#'              0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0), 16, 16, byrow = TRUE)
#' M <- c(3,5,6,15,16)
#' C <- c(4,7)
#' MAG(ex, M, C, plot=TRUE)
#' ###################################################
#' H <- matrix(c(0,100,1,0,100,0,100,0,0,100,0,100,0,1,100,0),4,4)
#' Max(H)
#'
`MAG`<-function (amat,M=c(),C=c(),showmat=TRUE,plot=FALSE, plotfun = plotGraph, ...)
{
	return(Max(AG(amat,M,C,showmat,plot, plotfun = plotGraph, ...)))
}
############################################################################
############################################################################
#Plot<-function(a)
#{
#	if(class(a)[1] == "igraph" || class(a)[1] == "graphNEL" || class(a)[1] == "character") {
#		a<-grMAT(a)}
#	if(is(a,"matrix")){
#		if(nrow(a)==ncol(a)){
#			if(length(rownames(a))!=ncol(a)){
#	 			rownames(a)<-1:ncol(a)
#	 			colnames(a)<-1:ncol(a)}
#			l1<-c()
#			l2<-c()
#			for (i in 1:nrow(a)){
#				for (j in i:nrow(a)){
#					if (a[i,j]==1){
#						l1<-c(l1,i,j)
#						l2<-c(l2,2)}
#					if (a[j,i]%%10==1){
#						l1<-c(l1,j,i)
#						l2<-c(l2,2)}
#					if (a[i,j]==10){
#						l1<-c(l1,i,j)
#						l2<-c(l2,0)}
#					if (a[i,j]==11){
#						l1<-c(l1,i,j,i,j)
#						l2<-c(l2,2,0)}
#					if (a[i,j]==100){
#						l1<-c(l1,i,j)
#						l2<-c(l2,3)}
#					if (a[i,j]==101){
#						l1<-c(l1,i,j,i,j)
#						l2<-c(l2,2,3)}
#					if (a[i,j]==110){
#						l1<-c(l1,i,j,i,j)
#						l2<-c(l2,0,3)}
#					if (a[i,j]==111){
#						l1<-c(l1,i,j,i,j,i,j)
#						l2<-c(l2,2,0,3)}
#					}
#				}
#			}
#		else {
#      	  stop("'object' is not in a valid adjacency matrix form")
#    	}
#    	if(length(l1)>0){
#		l1<-l1-1
#    		agr<-graph(l1,n=nrow(a),directed=TRUE)}
#    	if(length(l1)==0){
#    		agr<-graph.empty(n=nrow(a), directed=TRUE)
#		return(plot(agr,vertex.label=rownames(a)))}
#    	return( tkplot(agr, layout=layout.kamada.kawai, edge.curved=FALSE:TRUE, vertex.label=rownames(a),edge.arrow.mode=l2))
#	}
#	else {
#        stop("'object' is not in a valid format")}
#}
############################################################################
############################################################################


#' Markov equivalence for regression chain graphs.
#'
#' \code{MarkEqMag} determines whether two RCGs (or subclasses of RCGs) are
#' Markov equivalent.
#'
#' The function checks whether the two graphs have the same skeleton and
#' unshielded colliders.
#'
#' @param amat An adjacency matrix of an RCG or a graph that can be a
#' \code{graphNEL} or an \code{\link{igraph}} object or a vector of length
#' \eqn{3e}, where \eqn{e} is the number of edges of the graph, that is a
#' sequence of triples (type, node1label, node2label). The type of edge can be
#' \code{"a"} (arrows from node1 to node2), \code{"b"} (arcs), and \code{"l"}
#' (lines).
#' @param bmat The same as \code{amat}.
#' @return "Markov Equivalent" or "NOT Markov Equivalent".
#' @author Kayvan Sadeghi
#' @seealso \code{\link{MarkEqMag}}, \code{\link{msep}}
#' @references Wermuth, N. and Sadeghi, K. (2011). Sequences of regressions and
#' their independences. \emph{TEST}, To appear.
#' \url{http://arxiv.org/abs/1103.2523}.
#' @keywords graphs bidirected graph directed acyclic graph Markov equivalence
#' regression chain graph undirected graph multivariate
#' @examples
#'
#' H1<-matrix(c(0,100,0,0,0,100,0,100,0,0,0,100,0,0,0,1,0,0,0,100,0,0,1,100,0),5,5)
#' H2<-matrix(c(0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,100,0,0,1,100,0),5,5)
#' H3<-matrix(c(0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0),5,5)
#' #MarkEqRcg(H1,H2)
#' #MarkEqRcg(H1,H3)
#' #MarkEqRcg(H2,H3)
#'
MarkEqRcg<-function(amat,bmat)
	{
	if(class(amat)[1] == "igraph"){
		amat<-grMAT(amat)}
	if(class(amat)[1] == "graphNEL"){
		amat<-grMAT(amat)}
	if(class(amat)[1] == "character"){
		amat<-grMAT(amat)}
	if( length(rownames(amat))!=ncol(amat) | length(colnames(amat))!=ncol(amat)){
		rownames(amat)<-1:ncol(amat)
		colnames(amat)<-1:ncol(amat)}
	if(class(bmat)[1] == "igraph"){
		bmat<-grMAT(bmat)}
	if(class(bmat)[1] == "graphNEL"){
		bmat<-grMAT(bmat)}
	if(class(bmat)[1] == "character"){
		bmat<-grMAT(bmat)}
	if( length(rownames(bmat))!=ncol(bmat) | length(colnames(bmat))!=ncol(bmat)){
		rownames(bmat)<-1:ncol(bmat)
		colnames(bmat)<-1:ncol(bmat)}
	bmat<-bmat[rownames(amat),colnames(amat),drop=FALSE]
	na<-ncol(amat)
	nb<-ncol(bmat)
	if(na != nb){
		return(FALSE)}
	at<-which(amat+t(amat)+diag(na)==0,arr.ind=TRUE)
	bt<-which(bmat+t(bmat)+diag(na)==0,arr.ind=TRUE)
	if(identical(at,bt)==FALSE){
		return(FALSE)}
	ai<-c()
	bi<-c()
	if(dim(at)[1]!=0){
		for(i in 1:dim(at)[1]){
			for(j in 1:na){
				if((amat[at[i,1],j]%%10==1 || amat[at[i,1],j]>99) &&(amat[at[i,2],j]%%10==1 || amat[at[i,2],j]>99) ){
					ai<-c(ai,j)}
			if((bmat[bt[i,1],j]%%10==1 || bmat[bt[i,1],j]>99) &&(bmat[bt[i,2],j]%%10==1 || bmat[bt[i,2],j]>99) ){
				bi<-c(bi,j)}}
			if(identical(ai,bi)==FALSE){
				return(FALSE)}
		ai<-c()
		bi<-c()}}
	return(TRUE)
}
########################################################################################
########################################################################################


#' Markov equivalence of maximal ancestral graphs
#'
#' \code{MarkEqMag} determines whether two MAGs are Markov equivalent.
#'
#' The function checks whether the two graphs have the same skeleton and
#' colliders with order.
#'
#' @param amat An adjacency matrix of a MAG, or a graph that can be a
#' \code{graphNEL} or an \code{\link{igraph}} object or a vector of length
#' \eqn{3e}, where \eqn{e} is the number of edges of the graph, that is a
#' sequence of triples (type, node1label, node2label). The type of edge can be
#' \code{"a"} (arrows from node1 to node2), \code{"b"} (arcs), and \code{"l"}
#' (lines).
#' @param bmat The same as \code{amat}.
#' @return "Markov Equivalent" or "NOT Markov Equivalent".
#' @author Kayvan Sadeghi
#' @seealso \code{\link{MarkEqRcg}}, \code{\link{msep}}
#' @references Ali, R.A., Richardson, T.S. and Spirtes, P. (2009) Markov
#' equivalence for ancestral graphs. \emph{Annals of Statistics},
#' 37(5B),2808-2837.
#' @keywords graphs Markov equivalence maximal ancestral graphs multivariate
#' @examples
#'
#' H1<-matrix(  c(0,100,  0,  0,
#' 	         100,  0,100,  0,
#'                0,100,  0,100,
#'                0,  1,100,  0), 4, 4)
#' H2<-matrix(c(0,0,0,0,1,0,100,0,0,100,0,100,0,1,100,0),4,4)
#' H3<-matrix(c(0,0,0,0,1,0,0,0,0,1,0,100,0,1,100,0),4,4)
#' MarkEqMag(H1,H2)
#' MarkEqMag(H1,H3)
#' MarkEqMag(H2,H3)
#'
`MarkEqMag` <- function(amat,bmat)
	{
	if(class(amat)[1] == "igraph"){
		amat<-grMAT(amat)}
	if(class(amat)[1] == "graphNEL"){
		amat<-grMAT(amat)}
	if(class(amat)[1] == "character"){
		amat<-grMAT(amat)}
	if( length(rownames(amat))!=ncol(amat) | length(colnames(amat))!=ncol(amat)){
		rownames(amat)<-1:ncol(amat)
		colnames(amat)<-1:ncol(amat)}
	if(class(bmat)[1] == "igraph"){
		bmat<-grMAT(bmat)}
	if(class(bmat)[1] == "graphNEL"){
		bmat<-grMAT(bmat)}
	if(class(bmat)[1] == "character"){
		bmat<-grMAT(bmat)}
	if( length(rownames(bmat))!=ncol(bmat) | length(colnames(bmat))!=ncol(bmat)){
		rownames(bmat)<-1:ncol(bmat)
		colnames(bmat)<-1:ncol(bmat)}
	bmat<-bmat[rownames(amat),colnames(amat),drop=FALSE]
	na<-ncol(amat)
	nb<-ncol(bmat)
	if(na != nb){
		return(FALSE)}
	at<-which(amat+t(amat)+diag(na)==0,arr.ind=TRUE)
	bt<-which(bmat+t(bmat)+diag(na)==0,arr.ind=TRUE)
	if(identical(at,bt)==FALSE){
		return(FALSE)}
	a<-c()
	b<-c()
	if(length(at)>0){
		for(i in 1:dim(at)[1]){
			for(j in 1:na){
				if((amat[at[i,1],j]%%10==1 || amat[at[i,1],j]>99) &&(amat[at[i,2],j]%%10==1 || amat[at[i,2],j]>99) ){
					a<-rbind(a,c(at[i,1],j,at[i,2]))}
			if((bmat[bt[i,1],j]%%10==1 || bmat[bt[i,1],j]>99) &&(bmat[bt[i,2],j]%%10==1 || bmat[bt[i,2],j]>99) ){
					b<-rbind(b,c(bt[i,1],j,bt[i,2]))}}}
		if(identical(a,b)==FALSE){
			return(FALSE)}}
	ar<-which(amat%%10==1,arr.ind=TRUE)
	br<-which(bmat%%10==1,arr.ind=TRUE)
	ap<-c()
	bp<-c()
	if(length(ar)>0){
		for(i in 1:dim(ar)[1]){
			for(j in 1:na){
				if((amat[ar[i,1],j]>99) &&(amat[ar[i,2],j]>99)){
					ap<-rbind(ap,c(ar[i,1],j,ar[i,2]))}}}}
	if(length(br)>0){
		for(i in 1:dim(br)[1]){
			for(j in 1:nb){
				if((bmat[br[i,1],j]>99) &&(bmat[br[i,2],j]>99)){
					bp<-rbind(bp,c(br[i,1],j,br[i,2]))}}}}
	if(length(ap)>0){
		aptt<-ap
		apt<-c(ap,1)
		Qonen<-c()
		Qtwon<-c()
		while(length(apt)-length(aptt)>0){
			apt<-aptt
			for(i in (1:dim(ap)[1])){
				Qone<-ap[i,1]
				Qtwo<-ap[i,2]
				while(length(Qone)>0){
					for(l in (1:length(Qone))){
						J<-which(((amat+t(amat)+diag(na))[ap[i,3],]==0) & ((amat[Qone[l],]>99) | (amat[,Qone[l]]%%100==1)))
						for(j in J){
							for(k in 1:dim(a)[1]){
								if(min(a[k,]==c(j,Qone[l],Qtwo[l]))==1){
									a<-rbind(a,ap[i,])
									aptt<-aptt[(1:dim(aptt)[1])[-i],]
									break
									break
									break
									break}}}
						Q<-which((amat[,ap[i,3]]%%10==1) & (amat[Qone[l],]>99))
						for(q in Q){
							for(k in 1:dim(a)[1]){
								if(min(a[k,]==c(q,Qone[l],Qtwo[l]))==1){
									Qtwon<-c(Qtwon,Qone[l])
									Qonen<-c(Qonen,q)}}}}
					Qtwo<-Qtwon
					Qone<-Qonen
					Qonen<-c()
					Qtwon<-c()}}}}
	if(length(bp)>0){
		bptt<-bp
		bpt<-c(bp,1)
		Qbonen<-c()
		Qbtwon<-c()
		while(length(bpt)-length(bptt)>0){
			bpt<-bptt
			for(i in (1:dim(bp)[1])){
				Qbone<-bp[i,1]
				Qbtwo<-bp[i,2]
				while(length(Qbone)>0){
					for(l in (1:length(Qbone))){
						J<-which(((bmat+t(bmat)+diag(nb))[bp[i,3],]==0) & ((bmat[Qbone[l],]>99) | (bmat[,Qbone[l]]%%100==1)))
						for(j in J){
							for(k in 1:dim(b)[1]){
								if(min(b[k,]==c(j,Qbone[l],Qbtwo[l]))==1){
									b<-rbind(b,bp[i,])
									bptt<-bptt[(1:dim(bptt)[1])[-i],]
									break
									break
									break
									break}}}
						Qb<-which((bmat[,bp[i,3]]%%10==1) & (bmat[Qbone[l],]>99))
						for(q in Qb){
							for(k in 1:dim(b)[1]){
								if(min(b[k,]==c(q,Qbone[l],Qbtwo[l]))==1){
									Qbtwon<-c(Qbtwon,Qbone[l])
									Qbonen<-c(Qbonen,q)}}}}
					Qbtwo<-Qbtwon
					Qbone<-Qbonen
					Qbonen<-c()
					Qbtwon<-c()}}}}
	if(length(a)!=length(b)){
		return(FALSE)}
	f<-c()
	if((length(a)>0) && (length(b)>0)){
		for(i in 1:dim(a)[1]){
			for(j in 1:dim(b)[1]){
				f<-c(f,min(a[i,]==b[j,]))}
			if(max(f)==0){
				return(FALSE)}}}
	return(TRUE)
}
##########################################################################################
###########################################################################################


#' Representational Markov equivalence to undirected graphs.
#'
#' \code{RepMarUG} determines whether a given maximal ancestral graph can be
#' Markov equivalent to an undirected graph, and if that is the case, it finds
#' an undirected graph that is Markov equivalent to the given graph.
#'
#' \code{RepMarBG} looks for presence of an unshielded collider V-configuration
#' in graph.
#'
#' @param amat An adjacency matrix, or a graph that can be a \code{graphNEL} or
#' an \code{\link{igraph}} object or a vector of length \eqn{3e}, where \eqn{e}
#' is the number of edges of the graph, that is a sequence of triples (type,
#' node1label, node2label). The type of edge can be \code{"a"} (arrows from
#' node1 to node2), \code{"b"} (arcs), and \code{"l"} (lines).
#' @return A list with two components: \code{verify} and \code{amat}.
#' \code{verify} is a logical value, \code{TRUE} if there is a representational
#' Markov equivalence and \code{FALSE} otherwise.  \code{amat} is either
#' \code{NA} if \code{verify == FALSE} or the adjacency matrix of the generated
#' graph, if \code{verify == TRUE}. In this case it consists of 4 different
#' integers as an \eqn{ij}-element: 0 for a missing edge between \eqn{i} and
#' \eqn{j}, 1 for an arrow from \eqn{i} to \eqn{j}, 10 for a full line between
#' \eqn{i} and \eqn{j}, and 100 for a bi-directed arrow between \eqn{i} and
#' \eqn{j}. These numbers are added to be associated with multiple edges of
#' different types. The matrix is symmetric w.r.t full lines and bi-directed
#' arrows.
#' @author Kayvan Sadeghi
#' @seealso \code{\link{MarkEqMag}}, \code{\link{MarkEqRcg}},
#' \code{\link{RepMarBG}}, \code{\link{RepMarDAG}}
#' @references Sadeghi, K. (2011). Markov equivalences for subclasses of
#' loopless mixed graphs. \emph{Submitted}, 2011.
#' @keywords graphs bidirected graph Markov equivalence maximal ancestral graph
#' representational Markov equivalence
#' @examples
#'
#' H<-matrix(c(0,10,0,0,10,0,0,0,0,1,0,100,0,0,100,0),4,4)
#' RepMarUG(H)
#'
RepMarUG<-function(amat)
{
	if(class(amat)[1] == "igraph"){
		amat<-grMAT(amat)}
	if(class(amat)[1] == "graphNEL"){
		amat<-grMAT(amat)}
	if(class(amat)[1] == "character"){
		amat<-grMAT(amat)}
	if( length(rownames(amat))!=ncol(amat) | length(colnames(amat))!=ncol(amat)){
		rownames(amat)<-1:ncol(amat)
		colnames(amat)<-1:ncol(amat)}
	na<-ncol(amat)
	at<-which(amat+t(amat)+diag(na)==0,arr.ind=TRUE)
	if(dim(at)[1]!=0){
		for(i in 1:dim(at)[1]){
			for(j in 1:na){
				if((amat[at[i,1],j]%%10==1 || amat[at[i,1],j]>99) &&(amat[at[i,2],j]%%10==1 || amat[at[i,2],j]>99) ){
					return(list(verify = FALSE, amat = NA))}}}}
	for(i in 1:na){
		for(j in 1:na){
			if(amat[i,j]==100){
				amat[i,j]<-10}
			if(amat[i,j]==1){
				amat[i,j]<-10
				amat[j,i]<-10}}}
	return(list(verify = TRUE, amat = amat))
}
########################################################################################
#########################################################################################


#' Representational Markov equivalence to bidirected graphs.
#'
#' \code{RepMarBG} determines whether a given maximal ancestral graph can be
#' Markov equivalent to a bidirected graph, and if that is the case, it finds a
#' bidirected graph that is Markov equivalent to the given graph.
#'
#' \code{RepMarBG} looks for presence of an unshielded non-collider
#' V-configuration in graph.
#'
#' @param amat An adjacency matrix, or a graph that can be a \code{graphNEL} or
#' an \code{\link{igraph}} object or a vector of length \eqn{3e}, where \eqn{e}
#' is the number of edges of the graph, that is a sequence of triples (type,
#' node1label, node2label). The type of edge can be \code{"a"} (arrows from
#' node1 to node2), \code{"b"} (arcs), and \code{"l"} (lines).
#' @return A list with two components: \code{verify} and \code{amat}.
#' \code{verify} is a logical value, \code{TRUE} if there is a representational
#' Markov equivalence and \code{FALSE} otherwise.  \code{amat} is either
#' \code{NA} if \code{verify == FALSE} or the adjacency matrix of the generated
#' graph, if \code{verify == TRUE}. In this case it consists of 4 different
#' integers as an \eqn{ij}-element: 0 for a missing edge between \eqn{i} and
#' \eqn{j}, 1 for an arrow from \eqn{i} to \eqn{j}, 10 for a full line between
#' \eqn{i} and \eqn{j}, and 100 for a bi-directed arrow between \eqn{i} and
#' \eqn{j}. These numbers are added to be associated with multiple edges of
#' different types. The matrix is symmetric w.r.t full lines and bi-directed
#' arrows.
#' @author Kayvan Sadeghi
#' @seealso \code{\link{MarkEqMag}}, \code{\link{MarkEqRcg}},
#' \code{\link{RepMarDAG}}, \code{\link{RepMarUG}}
#' @references Sadeghi, K. (2011). Markov equivalences for subclasses of
#' loopless mixed graphs. \emph{Submitted}, 2011.
#' @keywords graphs bidirected graph Markov equivalence maximal ancestral graph
#' representational Markov equivalence
#' @examples
#'
#' H<-matrix(c(0,10,0,0,10,0,0,0,0,1,0,100,0,0,100,0),4,4)
#' RepMarBG(H)
#'
RepMarBG<-function(amat)
{
	if(class(amat)[1] == "igraph"){
		amat<-grMAT(amat)}
	if(class(amat)[1] == "graphNEL"){
		amat<-grMAT(amat)}
	if(class(amat)[1] == "character"){
		amat<-grMAT(amat)}
	if( length(rownames(amat))!=ncol(amat) | length(colnames(amat))!=ncol(amat)){
		rownames(amat)<-1:ncol(amat)
		colnames(amat)<-1:ncol(amat)}
	na<-ncol(amat)
	at<-which(amat+t(amat)+diag(na)==0,arr.ind=TRUE)
	if(dim(at)[1]!=0){
		for(i in 1:dim(at)[1]){
			for(j in 1:na){
				if(amat[at[i,1],j]%%100>9 || amat[j,at[i,1]]%%10==1 || amat[at[i,2],j]%%100>9 || amat[j,at[i,2]]%%10==1){
					return(list(verify= FALSE, amat = NA))}}}}
	for(i in 1:na){
		for(j in 1:na){
			if(amat[i,j]==10){
				amat[i,j]<-100}
			if(amat[i,j]==1){
				amat[i,j]<-100
				amat[j,i]<-100}}}
	return(list(verify = TRUE, amat = amat))
}
########################################################################################
########################################################################################


#' Representational Markov equivalence to directed acyclic graphs.
#'
#' \code{RepMarDAG} determines whether a given maximal ancestral graph can be
#' Markov equivalent to a directed acyclic graph, and if that is the case, it
#' finds a directed acyclic graph that is Markov equivalent to the given graph.
#'
#' \code{RepMarDAG} first looks whether the subgraph induced by full lines is
#' chordal and whether there is a minimal collider path or cycle of length 4 in
#' graph.
#'
#' @param amat An adjacency matrix, or a graph that can be a \code{graphNEL} or
#' an \code{\link{igraph}} object or a vector of length \eqn{3e}, where \eqn{e}
#' is the number of edges of the graph, that is a sequence of triples (type,
#' node1label, node2label). The type of edge can be \code{"a"} (arrows from
#' node1 to node2), \code{"b"} (arcs), and \code{"l"} (lines).
#' @return A list with two components: \code{verify} and \code{amat}.
#' \code{verify} is a logical value, \code{TRUE} if there is a representational
#' Markov equivalence and \code{FALSE} otherwise.  \code{amat} is either
#' \code{NA} if \code{verify == FALSE} or the adjacency matrix of the generated
#' graph, if \code{verify == TRUE}. In this case it consists of 4 different
#' integers as an \eqn{ij}-element: 0 for a missing edge between \eqn{i} and
#' \eqn{j}, 1 for an arrow from \eqn{i} to \eqn{j}, 10 for a full line between
#' \eqn{i} and \eqn{j}, and 100 for a bi-directed arrow between \eqn{i} and
#' \eqn{j}. These numbers are added to be associated with multiple edges of
#' different types. The matrix is symmetric w.r.t full lines and bi-directed
#' arrows.
#' @author Kayvan Sadeghi
#' @seealso \code{\link{MarkEqMag}}, \code{\link{MarkEqRcg}},
#' \code{\link{RepMarBG}}, \code{\link{RepMarUG}}
#' @references Sadeghi, K. (2011). Markov equivalences for subclasses of
#' loopless mixed graphs. \emph{Submitted}, 2011.
#' @keywords graphs bidirected graph Markov equivalence maximal ancestral graph
#' representational Markov equivalence
#' @examples
#'
#' H<-matrix(c(0,10,0,0,10,0,0,0,0,1,0,100,0,0,100,0),4,4)
#' RepMarBG(H)
#'
RepMarDAG<-function(amat)
{
	if(class(amat)[1] == "igraph"){
		amat<-grMAT(amat)}
	if(class(amat)[1] == "graphNEL"){
		amat<-grMAT(amat)}
	if(class(amat)[1] == "character"){
		amat<-grMAT(amat)}
	if(length(rownames(amat))!=ncol(amat) | length(colnames(amat))!=ncol(amat)){
		rownames(amat)<-1:ncol(amat)
		colnames(amat)<-1:ncol(amat)}
	na<-ncol(amat)
	full<-sort(unique(which(amat%%100>9,arr.ind=TRUE)[,1]))
	arc<-sort(unique(which(amat>99,arr.ind=TRUE)[,1]))
	arrow<-sort(unique(as.vector(which(amat%%10==1,arr.ind=TRUE))))
	S<-full[full!=full[1]]
	Ma<-full[1]
	while(length(S)>0){
		dim<-c()
		for(i in S){
			dim<-c(dim,length(which(amat[i,Ma]%%100>9,arr.ind=TRUE)))}
		s<-S[which(dim==max(dim))[1]]
		ns<-which(amat[s,Ma]%%100>9,arr.ind=TRUE)
		if(min(amat[ns,ns]+diag(length(ns)))==0){
			return(FALSE)}
		Ma<-c(Ma,s)
		S<-S[S!=s]}
	at<-which(amat+t(amat)==0,arr.ind=TRUE)
	ai<-c()
	if(length(at[,1])!=0){
		for(i in (1:length(at[,1]))){
			for(j in 1:na){
				if((amat[at[i,1],j]%%10==1 || amat[at[i,1],j]>99) &&(amat[at[i,2],j]%%10!=1) && (amat[at[i,2],j]<100) ){
					ai<-c(ai,j)}}
			if(length(ai)==0){
				break}
			for(j in ((1:na)[-ai])){
				if((max(amat[ai,j]>99)==1) &&(amat[at[i,2],j]%%10==1 || amat[at[i,2],j]>99) && (amat[at[i,1],j]%%10!=1) && (amat[at[i,1],j]<100)){
					return(list(verify = FALSE, amat = NA))}}
		ai<-c()}}
	for(i in Ma){
		v<-which(amat[i,]%%100>9)
		for(j in v){
			amat[i,j]<-1
			amat[j,i]<-0}}
	at<-which(amat+t(amat)+diag(na)==0,arr.ind=TRUE)
	if(length(at[,1])!=0){
		for(i in (1:length(at[,1]))){
			for(j in 1:na){
				if((amat[at[i,1],j]%%10==1 || amat[at[i,1],j]>99) &&(amat[at[i,2],j]%%10==1 || amat[at[i,2],j]>99) ){
					amat[at[i,1],j]<-1
					amat[at[i,2],j]<-1}}}}
	O<-c()
	Oc<-arrow
	while(length(Oc)>0){
		for(i in Oc){
	 		if(max(amat[i,Oc]%%10)==0){
		O<-c(O,i)
		break}}
		Oc<-Oc[Oc!=i]}
	for(i in arc){
		if(length(which(arrow==i))==0){
			O<-c(O,i)}}
	aarc<-which(amat>99,arr.ind=TRUE)
	if(length(aarc)[1]>0){
		for(i in 1:(length(aarc[,1]))){
		#	if(length(O)==0){
		#		amat[aarc[i,1],aarc[i,2]]<-0
			if(which(O==aarc[i,1])>which(O==aarc[i,2])){
				amat[aarc[i,1],aarc[i,2]]<-1}
			else{
				amat[aarc[i,1],aarc[i,2]]<-0}}}
	return(list(verify = TRUE, amat = amat))
}
##################################################################################
