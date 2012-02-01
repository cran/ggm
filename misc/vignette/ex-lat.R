### Example similar to that of Sewel et al. (1970)

rm(list = ls())

r = scan()
1     .288  .589  .438   .418
.288  1     .194  .359   .380
.589  .194  1     .473   .459 
.438  .359 .473    1     .611
.418  .380 .459   .611   1   


R = matrix(r, 5, 5)

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
library(ggm2)
#source("~/Documents/R/ggm/ggm2/R/sefit_senin.R")  

G = DAG(EA ~ AP+Z, SO ~ Z, Z ~ AP+MA+SES, SES ~MA, AP ~ MA)
out = fitDagLatent(G, R, 3500, latent = "Z", seed = 23, norm = 2)    
return()
sefit(G, "Z", out, R, 3500) 


## Simulation

source("~/Documents/R/ggm/ggm2/R/sefit_senin.R")  


a = scan()
1  0 -0.5 -0.05  0.00  0.0
0  1 -0.5  0.00  0.00  0.0
0  0  1.0  0.60  0.48  0.3
0  0  0.0  1.00  0.00  0.6
0  0  0.0  0.00  1.00  0.3
0  0  0.0  0.00  0.00  1.0

A = matrix(a, 6, 6, byrow = TRUE)
rownames(A) = c("EA",  "SO",    "Z" ,    "AP",    "SES",    "MA")
colnames(A) = c("EA",  "SO",    "Z" ,    "AP",    "SES",    "MA")  
                         
Delta = c(1,1,1,1,1,1)

X = rnormDag(3500, A, Delta) 
       
S = cov(X[,-3])   #remove latent
mo = fitDagLatent(G, S, 3500, norm= 2, latent = "Z", pri = 0) 
                   
sefit(G, "Z", mo, S, 3500)
                          
## Simulation of new data and ML fit

n = 100
B = rep(0, n)
for(i in 1:n){     
	y = rnormDag(n= 3500, A, Delta)  
	print(y[1:2,])        
    V = cov(y[,-3])   # remove latent        
	mo = fitDagLatent(G, V, 3500, norm= 2, latent = "Z", pri = 0)
    b =  mo$Ahat["Z", "AP"]
    B[i] = b
}

       sd(abs(B[abs(B) < 100]))

