########################### Examples ######################################################


ex<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ##The adjacency matrix of a DAG
0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,
1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,
0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
0,0,0,0,0,0,1,0,1,0,0),16,16,byrow=TRUE)

M<-c(3,5,6,15,16)
C<-c(4,7)
m<-c(M,13)
c<-c(C,8,14)
exgr<-as(ex,"graphNEL")   # not really needed 
exgr = ex

RG(exgr,c(),c()) #plotting the DAG     

RG(ex,M,C,plot=TRUE) #Plotting the ribbonless graph from the matrix
RG(ex,m,c) 
RG(RGMAT(ex,M,C),m,c) ##should be the same as the previous one 

SG(exgr,M,C) #ploting the summary graph from a graphNEL object
SG(ex,m,c)
SG(SGMAT(ex,M,C),8,c(3,9),byrownames=FALSE) ##this considers M and C subsets of {1,..,n} instead of the rownames, i.e the node labels, should be the same as the previous one 
AG(ex,M,C)
a<-c('b',1,2,'b',1,14,'a',9,8,'l',9,11,'a',10,8,'a',11,2,'a',11,10,'a',12,1,'b',12,14,'a',13,10,'a',13,12) ##another way of generating an IPG by defining the edge set
AG(a,c(),c()) ##should be the same a the previous one
AG(exgr,m,c)
msep(SGMAT(ex,M,C),13,c(9,14),8) ##Is 13 independent of {9,14} given 8?
msep(SGMAT(ex,M,C),13,c(9,14)) ##Is 13 independent of {9,14}?

AG(ex,m,c)
a<-c('b',1,2,'b',1,14,'a',9,8,'l',9,11,'a',10,8,'a',11,2,'a',11,10,'a',12,1,'b',12,14,'a',13,10,'a',13,12) 

rownames(ex)<-c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p")
colnames(ex)<-c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p")
M<-c("c","e","f","o","p")
C<-c("d","g")
RG(ex,M,C,plot=TRUE) #Plotting the ribbonless graph from the matrix

