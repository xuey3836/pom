library(MASS)
# library(rms)
library(maxLik)
library(Rcpp)

# setwd("/home/o0/Desktop/paper/0307/")
setwd("C:/Users/o0/Desktop/MSE for POM/code")
source("likelihood.R")

sourceCpp("likelihood.cpp")
sourceCpp("likelihoodgrad.cpp")
# sourceCpp("likelihoodmod.cpp")
# sourceCpp("likelihoodmodgrad.cpp")


#### functions defined:
#The function define in this R 

# phi                    1/(1+exp(-x))
# df.to.list             trasform a dataframe to list (need have name "y")
# p.likelihood           the pseudo-likelihood combine the lp and lm
# p.likelihood.grad      the gradient of  the  combine pseudo-likelihood
# p.likelihood.fisher    the fisher information of the  combine pseu-likelihood
# likelihood.mod         the modify likelihood 
# likelihood.mod.grad    the gradient of the modify likelihood
# likelihood             the pseudo-likelihood cobime the lp and lm using Rcpp
# likelihoodgrad         the gradient of  the  combine pseudo-likelihood using Rcpp
# cp.estimate            calculate the three estimtors and three p-value Using Rcpp
# p.estimate             calculate the three estimtors and three p-value Using R

####### The argument of cp.estimate() or p.estimate
# y: a vector with integer values representing levels (1 to J)
# G is a vector,contain 0,1,2
# x: a vector or matrix of numerical covariate values,the default is  null
# The lengths of y, G, x(or nrow(x) if matrix) must be the same.


#value return 

#for example
N=1000
y=sample(1:3,N,replace = T,prob = c(0.5,0.3,0.2))
G=sample(0:2,N,replace = T,prob = c(0.5,0.3,0.2))
x = rnorm(N,0,1)
# Sigma <- matrix(c(10,3,3,2),2,2)
# x=mvrnorm(n = N, rep(0, 2), Sigma)
system.time(cp.estimate(y=y,G=G,X=x))
system.time(p.estimate(y=y,G=G,X=x))
cp.estimate(y=y,G=G,X=x)
