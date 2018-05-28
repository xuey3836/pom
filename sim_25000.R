rm(list= ls())
library(MASS)
library(rms)
library(maxLik)
library(Rcpp)
library(netgwas)


setwd("/home/o0/Desktop/pom/code/")
# setwd("C:/Users/o0/Desktop/MSE for POM/code")
source("likelihood.R")

sourceCpp("likelihood.cpp")
sourceCpp("likelihoodgrad.cpp")

#--------------set parameter----------------
N = 100000
p = 100
n= c(500,300,200)
J = length(n)
alpha = c(3.48, 4.6)
beta = rep(0,p)
p0 =10
index0 = sample(1:p,p0,replace = FALSE)
# beta[index0]= runif(p*p0,-log(1.5),log(1.5))
beta[index0]= rep(log(1.4),p)
K = 1
gamma =rep(0.5,K)

#--------------generate the case control data for GWAS----------------
G=matrix(,nrow= N,ncol=p)
for (j in 1:ncol(G)){
  maf= runif(1,0.05,0.5)
  pG = c((1 - maf) ^ 2, 2 * maf * (1 - maf), maf ^ 2)
  G[,j] = sample(c(0,1,2),size = N,prob = pG,replace = T)
}
x = rnorm(N)
x = as.matrix(x)
py = (1 + exp(-outer(alpha, as.numeric(-G%*%beta - x%*%gamma)  , "+"))) ^ (-1)
aa = runif(N)
y = numeric(N)
for (i in 1:N)
  y[i] = sum(aa[i] > py[,i])
y = as.numeric(as.factor(y))
dat = data.frame(y = y, G = G, x = x)
#case control data
index1=sample(which(dat$y==1),n[1])
index2=sample(which(dat$y==2),n[2])
index3=sample(which(dat$y==3),n[3])
dat_case = dat[c(index1,index2,index3),]



#
y= sample(c(1:3),prob = c())


##--------------------simulation---------------------------------
ptm = proc.time()
sim = matrix(,nrow=3, ncol=p)
for (j in 1:p){
  if( j%%100== 0){
    cat("The loop is ", j,"\n")
  }
  y = dat_case$y
  x = dat_case$x
  g = dat_case[,j+1]
  sim[,j]=cp.estimate(y,G=g,X=x)[4:6]
}
proc.time() - ptm
##power
padj= apply(sim,1,p.adjust,method="BH")

