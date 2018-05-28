
library(MASS)
library(rms)
library(maxLik)
library(Rcpp)

# setwd("/home/o0/Desktop/paper/0307/")
setwd("C:/Users/o0/Desktop/MSE for POM/code")
source("likelihood.R")

sourceCpp("likelihood.cpp")
sourceCpp("likelihoodgrad.cpp")
sourceCpp("likelihoodmod.cpp")
sourceCpp("likelihoodmodgrad.cpp")

##---------------generate data function-----------------
genmaf.popdata = function(alpha, beta, gamma, maf, N) {
  pG = c((1 - maf) ^ 2, 2 * maf * (1 - maf), maf ^ 2)
  # x = NULL
  x=rnorm(N,0,1)
  # Sigma <- matrix(c(10,3,3,2),2,2)
  # x=mvrnorm(n = N, rep(0, 2), Sigma)
  
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
Pop.N = 30000000
#Pr(y=1)=0.97,P(Y=2)=0.02,P(Y=3)=0.01
##the number of category
J = 3
alpha = c(3.48, 4.6)
# alpha=c(4.18,5.27)
n = c(500, 300, 200)

# J=4
# alpha=c(3.89,4.59,5.51)
# n=c(1500,200,150,150)
##the dim of covariate
K = 1
gamma=rep(5,K)

##------------------------beta-------------------------
be = c(0, log(1.1), log(1.2), log(1.3), log(1.4), log(1.5))
# be=c(0,0.1,0.4,0.8,1,2)
ma = seq(0.05, 0.45, 0.05)

b.index = 1
beta = be[b.index]
#maf
m.index = 1
maf = ma[m.index]
##--------------------simulation---------------------------------
# ptm = proc.time()
# T = 1000
# Popdata=genmaf.popdata(alpha,beta,gamma,maf,N=Pop.N )
# Nc=unlist(lapply(Popdata,function(x){dim(x)[1]}))
# n=c(500,300,200)
# sim=replicate(T,{
#   cat("1")
#   ccdata = gen.ccdata(n, Popdata,Nc)
#   y = ccdata$y
#   G = ccdata$g
#   X = ccdata[,-(1:2)]
#   #c++
#   J = length(table(y))
#   N = length(y)
#   
#   K = ifelse(is.null(X), 0, ncol(as.matrix(X)))
#   #----------------- conbime likelihood-----------------
#   #--------------------------------pro and mod---------------------
#   pro.mod = lrm(y ~ ., data = ccdata)
#   p = pro.mod$coefficients
#   pro.var=pro.mod$var[J,J]
#   ui = matrix(0, 3 * J - 5, 3 * J + 2 * K - 1)
#   for (i in 1:(J - 2)) {
#     ui[i, i] = -1
#     ui[i, i + 1] = 1
#   }
#   ui[(J - 1):(2 * J - 3), (J + K + 1):(2 * J + K - 1)] = diag(J - 1)
#   for (i in (2 * J + K):(3 * J + K - 3)) {
#     ui[i - 2 - K, i] = -1
#     ui[i - 2 - K, i + 1] = 1
#   }
#   ci = rep(0, 3 * J - 5)
#   #init parameter
#   sta = c(3:(3 + J - 2), rep(1, K + 1), 5:(5 + J - 2), 1:(1 + J - 2), rep(1, K + 1))
#   ui.mod = ui[1:(2 * J - 3), 1:(2 * J + K - 1)]
#   ci.mod = ci[1:(2 * J - 3)]
#   sta.mod = sta[1:(2 * J + K - 1)]
#   mod.model = maxLik(
#     likelihoodmod,
#     grad = likelihoodmodgrad,
#     start = sta.mod,
#     constraints = list(ineqA = ui.mod, ineqB = ci.mod),
#     y = y,
#     X = as.matrix(ccdata[, -1]),
#     J = J,
#     lambda1=n[1]/N
#   )
#   m = mod.model$estimate
#   return(c(p,m))
#   #R
#   # p.estimate(y,G,X)
# })
# 
# proc.time() - ptm
# result=t(sim)
# value = result[, c(J,J+4)]
# ##bias
# bias=apply(value, 2, function(x) {mean(x)- beta})
# ##mse
# mse = apply(value, 2, function(x) {mean((x - beta) ^ 2)})
# 
# 
# boxplot(
#   result[, c(J,J+4)],
#   at = 1:2,
#   col = c(2,3),
#   boxwex = 0.5,
#   xlim = c(0, 4),
#   xaxt = "n",
#   xlab=expression(beta)
# )
# abline(h = beta, lty = 2)
# legend(3,-0.1, c("pro","mod"),pch=c(1,1,6),col=c(2,3),cex=1,text.font = 1)
# 
# boxplot(
#   result[,c(J+1,J+5)],
#   at = 1:2,
#   col = c(2,3),
#   boxwex = 0.5,
#   xlim = c(0, 4),
#   xaxt = "n",
#   xlab=expression(gamma)
# )
# abline(h = gamma, lty = 2)
# legend(3,4.5, c("pro","mod"),pch=c(1,1,6),col=c(2,3),cex=1,text.font = 1)

pm=NULL
T = 1000
for (i in 1:9){
  maf = ma[i] 
  cat("The maf is ", maf,"\n" )
  Popdata=genmaf.popdata(alpha,beta,gamma,maf,N=Pop.N )
  Nc=unlist(lapply(Popdata,function(x){dim(x)[1]}))
  n=c(500,300,200)
  sim=replicate(T,{
    ccdata = gen.ccdata(n, Popdata,Nc)
    y = ccdata$y
    G = ccdata$g
    X = ccdata[,-(1:2)]
    #c++
    J = length(table(y))
    N = length(y)
    
    K = ifelse(is.null(X), 0, ncol(as.matrix(X)))
    #----------------- conbime likelihood-----------------
    #--------------------------------pro and mod---------------------
    pro.mod = lrm(y ~ ., data = ccdata)
    p = pro.mod$coefficients
    pro.var=pro.mod$var[J,J]
    ui = matrix(0, 3 * J - 5, 3 * J + 2 * K - 1)
    for (i in 1:(J - 2)) {
      ui[i, i] = -1
      ui[i, i + 1] = 1
    }
    ui[(J - 1):(2 * J - 3), (J + K + 1):(2 * J + K - 1)] = diag(J - 1)
    for (i in (2 * J + K):(3 * J + K - 3)) {
      ui[i - 2 - K, i] = -1
      ui[i - 2 - K, i + 1] = 1
    }
    ci = rep(0, 3 * J - 5)
    #init parameter
    sta = c(3:(3 + J - 2), rep(1, K + 1), 5:(5 + J - 2), 1:(1 + J - 2), rep(1, K + 1))
    ui.mod = ui[1:(2 * J - 3), 1:(2 * J + K - 1)]
    ci.mod = ci[1:(2 * J - 3)]
    sta.mod = sta[1:(2 * J + K - 1)]
    mod.model = maxLik(
      likelihoodmod,
      grad = likelihoodmodgrad,
      start = sta.mod,
      constraints = list(ineqA = ui.mod, ineqB = ci.mod),
      y = y,
      X = as.matrix(ccdata[, -1]),
      J = J,
      lambda1=n[1]/N
    )
    m = mod.model$estimate
    return(c(p,m))
    #R
    # p.estimate(y,G,X)
  })
  result=t(sim)
  value = result[, c(J,J+4)]
  pm=cbind(value,pm)
  
}

write.csv(pm,"pm-no.csv")



