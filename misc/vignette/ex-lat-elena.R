### Example of simulation on Elena's DAG (2004)

rm(list = ls())



## Using fitDagLatent
require(ggm2)
source("~/Documents/R/ggm/ggm2/R/sefit_senin.R")  

G = DAG(EA ~ AP+Z, SO ~ Z, Z ~ AP+MA+SES, SES ~MA, AP ~ MA)
out = fitDagLatent(G, R, 3500, latent = "Z", seed = 23, norm = 2)
sefit(G, "Z", out, R, 3500) 


## Simulation for Elena's DAG

 
require(ggm)

a = scan()
     1  .4    0 -.2  0 
     0   1  -.8 1.2  0
     0   0    1   .6 0
     0   0    0    1  -.8
     0   0    0    0   1

A = matrix(a, 5, 5, byrow = TRUE) 

Delta = c(.8, .4, 1.1, 1, 0.9)


rownames(A) = c("X1",  "X2",    "X3" ,    "X4",    "X5")
colnames(A) =  c("X1",  "X2",    "X3" ,    "X4",    "X5")  
                         

X = data.frame(rnormDag(5000, A, Delta))    

  
S = cov(X[,-4])   #remove latent

G = DAG(X1 ~ X2+X4, X2 ~ X3 + X4, X3 ~ X4, X4 ~ X5)

# fitDagLatent takes a very long time to converge to dev = 0
# mo = fitDagLatent(G, S, 5000, norm= 2, latent = "X4", seed = 34) 
# source("~/Documents/R/ggm/ggm2/R/sefit_senin.R")  
# sefit(G, "X4", mo, S, 3500)
# The estimated se is wrong. 

             

# Correct estimate of a12
beta1 = lm(X1 ~ X2 + X3 + X5, data = X)$coef
beta2 = lm(X5 ~ X1 + X2 + X3, data = X)$coef  
a = beta1[2]; d = beta1[3];
e = beta2[2]; c = beta2[3]; b = beta2[4];

a12 = (a*b - c*d)/(b + d*e)

                   

## Simulation of new data and proper estimation of a12

samples = 1000 # number of simulations

A12 = rep(0, samples)
for(i in 1:samples){     
	X = data.frame(rnormDag(n= 5000, A, Delta))
	beta1 = lm(X1 ~ X2 + X3 + X5, data = X)$coef
	beta2 = lm(X5 ~ X1 + X2 + X3, data = X)$coef  
	a = beta1[2]; d = beta1[3];
	e = beta2[2]; c = beta2[3]; b = beta2[4];

	A12[i] = (a*b - c*d)/(b + d*e)

cat('.')
}   

#This is the right standard error

sd(A12)
   
## Now try with sem

require(sem)


mod = specify.model()
X1 <- X2, 1.2, NA
X1 <- X4, 1.4, NA
X2 <- X3, 2.3, NA
X2 <- X4, 2.4, NA
X3 <- X4, 3.4, NA
X4 <- X5, 4.5 , NA
X1 <-> X1, 1.1, NA
X2 <-> X2, 2.2, NA
X3 <-> X3, 3.3, NA
X4 <-> X4, NA, 1
X5 <-> X5, 5.5, NA


summary(sem(mod, S, 5000))

