### Example similar to that of Sewel et al. (1970)

rm(list = ls())

R = matrix(c(
1    , .288,  .589,  .438,   .418,
.288 , 1   ,  .194,  .359,   .380,
.589 , .194,  1   ,  .473,   .459, 
.438 , .359, .473 ,   1  ,   .611,
.418 , .380, .459 ,  .611,   1), 5, 5)   




require(sem)


nam = c("MA", "SES", "AP", "SO", "EA")
dimnames(R) = list(nam, nam)

mod = specifyModel()
EA <- AP, EA.AP, NA
EA <- Z, EA.Z, NA
SO <- Z, SO.Z, NA
Z <- AP, Z.AP, NA
Z <- SES, Z.SES, NA
Z <- MA, Z.MA,  1
AP <- MA, AP.MA, NA
SES <- MA, SES.MA, NA
EA <-> EA, EA.EA, NA
AP <-> AP, AP.AP, NA
Z<-> Z, NA, 1
MA <-> MA, MA.MA, 1
SES <-> SES, SES.SES, NA
SO <-> SO, SO.SO, NA




summary(sem(mod, R, 3500))


## Using fitDagLatent
require(ggm2)


G = DAG(EA ~ AP+Z, SO ~ Z, Z ~ AP+MA+SES, SES ~MA, AP ~ MA)
out = fitDagLatent(G, R, 3500, latent = "Z", norm = 2)    


`expa` <- function (B) 
{
### Finds the matrix P such that vec(B) = P b     	
	notzero <- as.vector(B) != 0
	d <- sum(notzero)	
	P <- matrix(0, prod(dim(B)), d)
 	P[notzero,] <- diag(d)
    P 
}
               
A = out$Ahat
B = diag(6) - A
P = expa(B)
D = diag(out$Dhat)     
Q = expa(D); Q = Q[,-6]
S = out$Shat
Ibb = t(P) %*% (S %x% solve(D)) %*% P
Ibw = t(P) %*% (solve(A) %x% solve(D)) %*% Q
Iww = 0.5 * t(Q) %*% (solve(D) %x% solve(D)) %*% Q
I = rbind(cbind(Ibb, Ibw), cbind(t(Ibw), Iww))
V = solve(I) / 3500
SE = sqrt(diag(V))

`allEdges` <-   function(amat, sep = "<-"){
### Finds all the edges of a graph with edge matrix amat.
    nn <- rownames(amat)
    E <- c()
    if(all(amat == t(amat))) { 
      amat[lower.tri(amat)] <- 0
    }
    for(i in nn) {
      e <- nn[amat[i,] == 1]
      if(length(e) == 0) next
      li <- paste(i, e, sep = sep)
      #dimnames(li) <- list(rep("", length(e)), rep("", 2))
      E <- c(E, li) 
    }
    E
  }
n1 = allEdges(In(B), sep = "<-")
n2 = allEdges(S*0 + diag(6), sep = "<->")
names(SE) = c(n1, n2[-6])
SE

