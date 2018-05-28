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
# sourceCpp("likelihoodmod.cpp")
# sourceCpp("likelihoodmodgrad.cpp")

##---------------generate data function-----------------
genmaf.popdata = function(alpha, beta, gamma, maf, N) {
  pG = c((1 - maf) ^ 2, 2 * maf * (1 - maf), maf ^ 2)
  k= length(gamma)
  if (k==0){
    x=NULL
  }else if(k==1){
    x=rnorm(N,0,1)
  }else if(k==2){
    Sigma <- matrix(c(10,3,3,2),2,2)
    x=mvrnorm(n = N, rep(0, 2), Sigma)
  }
  
  g = sample(0:2, N, replace = TRUE, prob = pG)
  y = numeric(N)
  
  if (is.null(x)) {
    py = (1 + exp(-outer(alpha, -beta * g  , "+"))) ^ (-1)
    aa = runif(N)
    for (i in 1:N)
      y[i] = sum(aa[i] > py[, i])
    y = as.numeric(as.factor(y))
    data=data.frame(y = y, g = g)
  }
  else{
    x = as.matrix(x)
    py = (1 + exp(-outer(alpha, -beta * g - rowSums(gamma * x)  , "+"))) ^(-1)
    aa = runif(N)
    for (i in 1:N)
      y[i] = sum(aa[i] > py[, i])
    y = as.numeric(as.factor(y))
    data=data.frame(y = y, g = g, x = x)
  }
  pop=list()
  for (j in 1:J){
    pop[[j]] = subset(data, y == j)
  }
  return(pop)
}

gen.ccdata = function(n, data,Nc) {
  cs <- NULL
  for (j in 1:J) {
    new <- data[[j]][sample(1:Nc[j], n[j]), ]
    cs = rbind(cs, new)
  }
  return(cs)
}


#------------set the parameter------------------------
Pop.N = 300000
#Pr(y=1)=0.97,P(Y=2)=0.02,P(Y=3)=0.01
##the number of category
J = 3
alpha = c(3.48, 4.6)
# alpha=c(4.18,5.27)


# J=4
# alpha=c(3.89,4.59,5.51)
# n=c(1500,200,150,150)
##the dim of covariate
K = 1
gamma=rep(0.5,K)

##------------------------beta-------------------------
be = c(0, log(1.1), log(1.2), log(1.3), log(1.4), log(1.5))
# be=c(0,0.1,0.4,0.8,1,2)

ma = seq(0.05, 0.45, 0.05)

b.index = 5
beta = be[b.index]
#maf
m.index =7
maf = ma[m.index]
##--------------------simulation---------------------------------
ptm = proc.time()
T = 25000
n= c(500,300,200)*3
sim = matrix(,T,3)
for (i in 1:T){
  if(i%%1000 == 0){cat("the loop is ", i,"\n")}
  Popdata=genmaf.popdata(alpha,beta,gamma,maf,N=Pop.N )
  Nc=unlist(lapply(Popdata,function(x){dim(x)[1]}))
  ccdata = gen.ccdata(n, Popdata,Nc)
  y = ccdata$y
  G = ccdata$g
  X = ccdata[,-(1:2)]
  #c++
  sim[i,]=cp.estimate(y,G,X)[4:6] 
}
proc.time()-ptm

##power
por =apply(sim < 0.05/T, 2, mean, rm.na = TRUE)
maf
n
beta
por
write.csv(por,paste0("por",m.index,b.index,".csv"))


