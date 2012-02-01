# Perugia, 11/08/2011
# An R routine to calculate the standard error using 
# missing information principle (EM Algorithm)
# Last modified: 11/08/2011 08:08:46
#

sefit <- function(amat, latent, fit, Syy, n)
{
  cmqi <- function (Syy, Sigma, z)
        {                                          
        y <- !z
        Q <- solve(Sigma)
        Qzz <- Q[z,z]
        Qzy <- Q[z,y]
        B   <- - solve(Qzz) %*% Qzy
        BSyy <- B %*% Syy
        E <- Sigma*0
        E[y,y] <- Syy
        E[y,z] <- t(BSyy)
        E[z,y] <- BSyy
        E[z,z] <- BSyy %*% t(B) + solve(Qzz)
        dimnames(E) <- dimnames(Sigma)
        E
        } 
A     <-(fit$Ahat)                             
Shat  <- (fit$Shat)                            
Khat  <- solve(Shat)                          
Dhat  <- (fit$Dhat)                            
Delta <- diag(fit$Dhat)
d1    <-solve(Delta)
dimnames(d1)   <- dimnames(A)
dimnames(Delta)<-dimnames(d1)
AA   <- solve(A)
nod  <- rownames(amat)
nam  <- rownames(Syy)
sek  <- intersect(nam, nod)
sek  <- c(sek, latent)
amat <- amat[sek,sek, drop=FALSE]
nod<-rownames(amat)
wherez <- is.element(nod, latent)
QQ <-cmqi(Syy, Shat, wherez)
q <- ncol(A)
H <- matrix(0, q ,q)
H[wherez,wherez] <- 1/(Khat[wherez,wherez])
Qtil  <- QQ-H
e21_1 <- Qtil %*% t(A) %*% d1     
e21_2 <- QQ %*% t(A) %*% d1       
e24_1 <- A %*% QQ %*% t(A)
e24_2 <- A %*% Qtil %*% t(A)
ij <- matrix(nod[allEdges(amat)], ncol=2)
ij[, 2:1]
p <- nrow(ij)
k <- c()
#----------------------------------------------------------------------------
for(u in 1:p)
  {
    for(v in 1:p)
      {
        i <- ij[u,1]; j <- ij[u,2]
        l <- ij[v,1]; m <- ij[v,2]
        new1 <-   QQ[j,m]%*%d1[l,i] + AA[m,i]%*%AA[j,l] + 2%*%e21_1[m,i]%*% e21_1[j,l] - 2%*%e21_2[m,i]%*% e21_2[j,l]
        k = c(k,new1)
        }
    }               
#----------------------------------------------------------------------------
k <- matrix(k, p, p, byrow=TRUE)
ed <- paste(ij[,1],"->", ij[,2], sep="")
dimnames(k) <- list(ed,ed)
p <- ncol(Syy)
p1 <-ncol(Dhat)
npar <- apply(topSort(amat), 2, sum)
kinv <- solve(k)
seA <- sqrt(diag(kinv))
# for delta
uu= 0.5*(Delta %*%Delta) - 0.5* (e24_1%*% e24_1) + (1/2) *( e24_2%*% e24_2 )
seD=diag(sqrt((uu)))
DD=Dhat/seD
DD=as.vector(DD)
P=2*(1-pnorm(abs(DD)))
tab2 = cbind(Estimated=Dhat,s.e.=seD,z_val=DD,p_val=P)
cat("\n")
n<-ncol(A)
pp = c()
At<- -t(A)
for (i in 1:n) 
  {
      for (j in 1:n)  
        {
          if (amat[i,j]==1)
              pp=cbind(pp,At[i,j])
        }
  } 
  
pp <-as.vector(pp)
TT <- pp/seA
TT <- as.vector(TT)
names(TT) <- names(seA)                                       
PP<-2*(1-pnorm(abs(TT)))
tab <-round(cbind(Estimated = pp, s.e. = seA, z_val = TT, p_val = PP),6)                               
list(seA=tab,seD=tab2)
}
