#The function define in this R 

# phi                    1/(1+exp(-x))
# df.to.list             trasform a dataframe to list (need have name "y")
# p.likelihood           the pseu-likelihood cobime the lp and lm
# p.likelihood.grad      the gradient of  the  combine pseu-likelihood
# p.likelihood.fisher    the fisher information of the  combine pseu-likelihood
# cp.estimate            estimate the three estimtors and three p-value Using Rcpp
# p.estimate             estimate the three estimtors and three p-value Using R


#y=1/(1+exp(-x))
phi <- function(x) {
  return((1 + exp(-x)) ^ (-1))
}
# trasform dataframe ,must have "y",to a list by ordinal "y"
df.to.list = function(df_data) {
  J = length(table(df_data$y))
  cs = list()
  for (j in 1:J) {
    cs[[j]] = as.matrix(subset(df_data, y == j)[, -1])
  }
  return(cs)
}
# the likelihood function for combine lp and lm
p.likelihood = function(theta, data) {
  J = length(data)
  
  K = dim(as.matrix(data[[1]]))[2] - 1
  
  beta <- theta[J:(J + K)]
  
  alpha <- c(-Inf, theta[1:(J - 1)], Inf)
  
  lambda0 <- theta[(J + K + 1):(2 * J + K - 1)]
  
  alpha.p <- c(-Inf, theta[(2 * J + K):(3 * J + K - 2)], Inf)
  beta.p <- theta[-(1:(3 * J + K - 2))]
  
  n <- unlist(lapply(data, function(x) {
    dim(x)[1]
  }))
  Lp = 0
  Lm = 0
  N = sum(n)
  lambda <- c(n[1] / N, lambda0)
  for (j in 1:J) {
    Lp = Lp + sum(log(phi(alpha.p[j + 1] - beta.p %*% t(data[[j]])) - phi(alpha.p[j] -
                                                                            beta.p %*% t(data[[j]]))))
  }
  for (j in 1:J) {
    a = sum(log(lambda[j] * (
      phi(alpha[j + 1] - beta %*% t(data[[j]])) - phi(alpha[j] - beta %*% t(data[[j]]))
    )))
    b = apply(data[[j]], 1, function(x) {
      log(sum(lambda * (
        phi(alpha - beta %*% x)[2:(J + 1)] - phi(alpha - beta %*% x)[1:J]
      )))
    })
    Lm = Lm + a - sum(b)
  }
  return(Lp + Lm)
}
##the first gradient of the likelihood function for combine lp and lm
p.likelihood.grad = function(theta, data) {
  J = length(data)
  
  K = dim(as.matrix(data[[1]]))[2] - 1
  
  beta <- theta[J:(J + K)]
  
  alpha <- c(-Inf, theta[1:(J - 1)], Inf)
  
  lambda0 <- theta[(J + K + 1):(2 * J + K - 1)]
  
  alpha.p <- c(-Inf, theta[(2 * J + K):(3 * J + K - 2)], Inf)
  beta.p <- theta[-(1:(3 * J + K - 2))]
  
  n <- unlist(lapply(data, function(x) {
    dim(x)[1]
  }))
  
  N = sum(n)
  lambda <- c(n[1] / N, lambda0)
  
  galn = rep(0, K + 1)
  gbln = rep(0, J - 1)
  gcln = rep(0, J - 1)
  #beta
  
  for (k in 1:(K + 1)) {
    for (j in 1:J) {
      galn[k] = galn[k] + sum(-(1 - phi(alpha[j + 1] - beta %*% t(data[[j]])) -
                                  phi(alpha[j] - beta %*% t(data[[j]]))) * data[[j]][, k])
      for (i in 1:n[j]) {
        x = data[[j]][i, ]
        a1 = sum(lambda * (1 - phi(alpha - beta %*% x)[2:(J + 1)] - phi(alpha -
                                                                          beta %*% x)[1:J]) * (phi(alpha - beta %*% x)[2:(J + 1)] - phi(alpha - beta %*%
                                                                                                                                          x)[1:J]))
        b1 = sum(lambda * (phi(alpha - beta %*% x)[2:(J + 1)] - phi(alpha -
                                                                      beta %*% x)[1:J]))
        galn[k] = galn[k] + x[k] * a1 / b1
      }
    }
  }
  
  for (j in (1:(J - 1))) {
    a2 = sum(phi(alpha[j + 1] - beta %*% t(data[[j]])) * (1 - phi(alpha[j +
                                                                          1] - beta %*% t(data[[j]]))) / (phi(alpha[j + 1] - beta %*% t(data[[j]])) -
                                                                                                            phi(alpha[j] - beta %*% t(data[[j]]))))
    b2 = sum(phi(alpha[j + 1] - beta %*% t(data[[j + 1]])) * (1 - phi(alpha[j +
                                                                              1] - beta %*% t(data[[j + 1]]))) / (phi(alpha[j + 2] - beta %*% t(data[[j +
                                                                                                                                                        1]])) - phi(alpha[j + 1] - beta %*% t(data[[j + 1]]))))
    gbln[j] = a2 - b2
    gcln[j] = n[j + 1] / lambda[j + 1]
    
    for (l in (1:J)) {
      for (i in 1:n[l]) {
        x = data[[l]][i, ]
        a3 = phi(alpha[j + 2] - beta %*% x) - phi(alpha[j + 1] - beta %*%
                                                    x)
        b3 = (lambda[j] - lambda[j + 1]) * phi(alpha[j + 1] - beta %*% x) *
          (1 - phi(alpha[j + 1] - beta %*% x))
        c3 = sum(lambda * (phi(alpha - beta %*% x)[2:(J + 1)] - phi(alpha -
                                                                      beta %*% x)[1:J]))
        gbln[j] = gbln[j] - b3 / c3
        gcln[j] = gcln[j] - a3 / c3
      }
    }
  }
  
  ##p
  galn.p = rep(0, K + 1)
  gbln.p = rep(0, J - 1)
  for (k in 1:(K + 1)) {
    for (j in 1:J) {
      galn.p[k] = galn.p[k] + sum(-(
        1 - phi(alpha.p[j + 1] - beta.p %*% t(data[[j]])) - phi(alpha.p[j] - beta.p %*%
                                                                  t(data[[j]]))
      ) * data[[j]][, k])
    }
  }
  
  for (j in (1:(J - 1))) {
    a2 = sum(phi(alpha.p[j + 1] - beta.p %*% t(data[[j]])) * (1 - phi(alpha.p[j +
                                                                                1] - beta.p %*% t(data[[j]]))) / (phi(alpha.p[j + 1] - beta.p %*% t(data[[j]])) -
                                                                                                                    phi(alpha.p[j] - beta.p %*% t(data[[j]]))))
    b2 = sum(phi(alpha.p[j + 1] - beta.p %*% t(data[[j + 1]])) * (1 - phi(alpha.p[j +
                                                                                    1] - beta.p %*% t(data[[j + 1]]))) / (phi(alpha.p[j + 2] - beta.p %*% t(data[[j +
                                                                                                                                                                    1]])) - phi(alpha.p[j + 1] - beta.p %*% t(data[[j + 1]]))))
    gbln.p[j] = a2 - b2
    
  }
  
  return(c(gbln, galn, gcln, gbln.p, galn.p))
}

# the fisher information matrix for the combined likelihood 
p.likelihood.fisher = function(theta, data) {
  J = length(data)
  
  K = dim(as.matrix(data[[1]]))[2] - 1
  
  beta <- theta[J:(J + K)]
  
  alpha <- c(-Inf, theta[1:(J - 1)], Inf)
  
  lambda0 <- theta[(J + K + 1):(2 * J + K - 1)]
  
  alpha.p <- c(-Inf, theta[(2 * J + K):(3 * J + K - 2)], Inf)
  beta.p <- theta[-(1:(3 * J + K - 2))]
  
  n <- unlist(lapply(data, function(x) {
    dim(x)[1]
  }))
  N = sum(n)
  lambda <- c(n[1] / N, lambda0)
  
  
  I = matrix(0, 3 * J - 1 + 2 * K, 3 * J - 1 + 2 * K)
  
  for (j in 1:J) {
    for (i in 1:n[j]) {
      x = data[[j]][i, ]
      ph.j = phi(alpha[j + 1] - beta %*% x)
      ph.jm = phi(alpha[j] - beta %*% x)
      
      ##pro
      ph.j.p = phi(alpha.p[j + 1] - beta.p %*% x)
      ph.jm.p = phi(alpha.p[j] - beta.p %*% x)
      #beta
      beta.grad = -x * (1 - ph.j - ph.jm)
      a = x * sum(lambda * (1 - phi(alpha[2:(J + 1)] - beta %*% x) - phi(alpha[1:J] -
                                                                           beta %*% x)) * (phi(alpha[2:(J + 1)] - beta %*% x) - phi(alpha[1:J] - beta %*%
                                                                                                                                      x)))
      b = sum(lambda * (phi(alpha[2:(J + 1)] - beta %*% x) - phi(alpha[1:J] -
                                                                   beta %*% x)))
      beta.grad = beta.grad + a / b
      
      ##pro
      beta.grad.p = -x * (1 - ph.j.p - ph.jm.p)
      #theta
      alpha.grad = rep(0, J - 1)
      lambda.grad = rep(0, J)
      
      alpha.grad.p = rep(0, J - 1)
      if (j == 1) {
        alpha.grad[j] = ph.j * (1 - ph.j) / (ph.j - ph.jm) - (lambda[j] - lambda[j +
                                                                                   1]) * ph.j * (1 - ph.j) / b
        alpha.grad.p[j] = ph.j.p * (1 - ph.j.p) / (ph.j.p - ph.jm.p)
        for (t in 2:(J - 1)) {
          ph.t = phi(alpha[t + 1] - beta %*% x)
          alpha.grad[t] = -(lambda[t] - lambda[t + 1]) * ph.t * (1 - ph.t) /
            b
        }
        for (t in 2:J) {
          ph.t = phi(alpha[t + 1] - beta %*% x)
          ph.tm = phi(alpha[t] - beta %*% x)
          lambda.grad[t] = -(ph.t - ph.tm) / b
        }
        #p
        
        theta.grad = c(alpha.grad,
                       beta.grad,
                       lambda.grad[-1],
                       alpha.grad.p,
                       beta.grad.p)
        I = I + theta.grad %*% t(theta.grad)
      }
      else{
        if (j < J) {
          alpha.grad[j] = ph.j * (1 - ph.j) / (ph.j - ph.jm) - (lambda[j] - lambda[j +
                                                                                     1]) * ph.j * (1 - ph.j) / b
          alpha.grad[j - 1] = -ph.jm * (1 - ph.jm) / (ph.j - ph.jm) - (lambda[j -
                                                                                1] - lambda[j]) * ph.jm * (1 - ph.jm) / b
          
          alpha.grad.p[j] = ph.j.p * (1 - ph.j.p) / (ph.j.p - ph.jm.p)
          alpha.grad.p[j - 1] = -ph.jm.p * (1 - ph.jm.p) / (ph.j.p - ph.jm.p)
          
          s.alpha = setdiff(c(1:(J - 1)), c(j - 1, j))
          for (t in s.alpha) {
            ph.t = phi(alpha[t + 1] - beta %*% x)
            alpha.grad[t] = -(lambda[t] - lambda[t + 1]) * ph.t * (1 - ph.t) /
              b
          }
          lambda.grad[j] = 1 / lambda[j] - (ph.j - ph.jm) / b
          s.lambda = setdiff(2:J, j)
          for (t in s.lambda) {
            ph.t = phi(alpha[t + 1] - beta %*% x)
            ph.tm = phi(alpha[t] - beta %*% x)
            lambda.grad[t] = -(ph.t - ph.tm) / b
          }
          theta.grad = c(alpha.grad,
                         beta.grad,
                         lambda.grad[-1],
                         alpha.grad.p,
                         beta.grad.p)
          I = I + theta.grad %*% t(theta.grad)
        }
        else{
          alpha.grad[j - 1] = -ph.jm * (1 - ph.jm) / (ph.j - ph.jm) - (lambda[j -
                                                                                1] - lambda[j]) * ph.jm * (1 - ph.jm) / b
          alpha.grad.p[j - 1] = -ph.jm.p * (1 - ph.jm.p) / (ph.j.p - ph.jm.p)
          for (t in 1:(J - 2)) {
            ph.t = phi(alpha[t + 1] - beta %*% x)
            alpha.grad[t] = -(lambda[t] - lambda[t + 1]) * ph.t * (1 - ph.t) /
              b
          }
          lambda.grad[j] = 1 / lambda[j] - (ph.j - ph.jm) / b
          for (t in 2:(J - 1)) {
            ph.t = phi(alpha[t + 1] - beta %*% x)
            ph.tm = phi(alpha[t] - beta %*% x)
            lambda.grad[t] = -(ph.t - ph.tm) / b
          }
          theta.grad = c(alpha.grad,
                         beta.grad,
                         lambda.grad[-1],
                         alpha.grad.p,
                         beta.grad.p)
          I = I + theta.grad %*% t(theta.grad)
        }
      }
    }
  }
  return(I / N)
}
likelihood.mod = function(theta, data) {
  J = length(data)
  
  K = dim(as.matrix(data[[1]]))[2] - 1
  
  beta <- theta[J:(J + K)]
  
  alpha <- c(-Inf, theta[1:(J - 1)], Inf)
  
  lambda0 <- theta[(J + K + 1):(2 * J + K - 1)]
  n <- unlist(lapply(data, function(x) {
    dim(x)[1]
  }))
  Lm = 0
  N = sum(n)
  lambda <- c(n[1] / N, lambda0)
  
  for (j in 1:J) {
    a = sum(log(lambda[j] * (
      phi(alpha[j + 1] - beta %*% t(data[[j]])) - phi(alpha[j] - beta %*% t(data[[j]]))
    )))
    b = apply(data[[j]], 1, function(x) {
      log(sum(lambda * (
        phi(alpha - beta %*% x)[2:(J + 1)] - phi(alpha - beta %*% x)[1:J]
      )))
    })
    Lm = Lm + a - sum(b)
  }
  return(Lm)
}
# lm=function(theta,y,X){
#   n=table(y)
#   J=length(n)
#   K=dim(X)[2]
#
#   beta <- theta[J:(J+K-1)]
#
#   alpha<- c(-Inf, theta[1:(J-1)], Inf)
#
#   lambda0 <- theta[(J+K):(2*J+K-2)]
#
#   Lm=0
#   N=sum(n)
#   lambda <- c(n[1]/N,lambda0)
#   for (i in 1:N){
#     betax=beta%*%X[i,]
#     for(j in 1:J){
#       if(y[i]==j){
#        a=lambda[j]*(1/(1+exp(-alpha[j+1]+betax)) - 1/(1+exp(-alpha[j]+betax)) );
#        b=0;
#         for (k in 1:J)
#         {
#           b=b+lambda[k]*( 1/(1+exp(-alpha[k+1]+betax)) - 1/(1+exp(-alpha[k]+betax)) );
#         }
#         Lm=Lm+log(a/b);
#       }
#       }
#   }
#   return(Lm)
# }

likelihood.mod.grad = function(theta, data) {
  J = length(data)
  
  K = dim(as.matrix(data[[1]]))[2] - 1
  
  beta <- theta[J:(J + K)]
  
  alpha <- c(-Inf, theta[1:(J - 1)], Inf)
  
  lambda0 <- theta[(J + K + 1):(2 * J + K - 1)]
  
  
  n <- unlist(lapply(data, function(x) {
    dim(x)[1]
  }))
  
  N = sum(n)
  lambda <- c(n[1] / N, lambda0)
  
  galn = rep(0, K + 1)
  gbln = rep(0, J - 1)
  gcln = rep(0, J - 1)
  #beta
  
  for (k in 1:(K + 1)) {
    for (j in 1:J) {
      galn[k] = galn[k] + sum(-(1 - phi(alpha[j + 1] - beta %*% t(data[[j]])) -
                                  phi(alpha[j] - beta %*% t(data[[j]]))) * data[[j]][, k])
      for (i in 1:n[j]) {
        x = data[[j]][i, ]
        a1 = sum(lambda * (1 - phi(alpha - beta %*% x)[2:(J + 1)] - phi(alpha -
                                                                          beta %*% x)[1:J]) * (phi(alpha - beta %*% x)[2:(J + 1)] - phi(alpha - beta %*%
                                                                                                                                          x)[1:J]))
        b1 = sum(lambda * (phi(alpha - beta %*% x)[2:(J + 1)] - phi(alpha -
                                                                      beta %*% x)[1:J]))
        galn[k] = galn[k] + x[k] * a1 / b1
      }
    }
  }
  
  for (j in (1:(J - 1))) {
    a2 = sum(phi(alpha[j + 1] - beta %*% t(data[[j]])) * (1 - phi(alpha[j +
                                                                          1] - beta %*% t(data[[j]]))) / (phi(alpha[j + 1] - beta %*% t(data[[j]])) -
                                                                                                            phi(alpha[j] - beta %*% t(data[[j]]))))
    b2 = sum(phi(alpha[j + 1] - beta %*% t(data[[j + 1]])) * (1 - phi(alpha[j +
                                                                              1] - beta %*% t(data[[j + 1]]))) / (phi(alpha[j + 2] - beta %*% t(data[[j +
                                                                                                                                                        1]])) - phi(alpha[j + 1] - beta %*% t(data[[j + 1]]))))
    gbln[j] = a2 - b2
    gcln[j] = n[j + 1] / lambda[j + 1]
    
    for (l in (1:J)) {
      for (i in 1:n[l]) {
        x = data[[l]][i, ]
        a3 = phi(alpha[j + 2] - beta %*% x) - phi(alpha[j + 1] - beta %*%
                                                    x)
        b3 = (lambda[j] - lambda[j + 1]) * phi(alpha[j + 1] - beta %*% x) *
          (1 - phi(alpha[j + 1] - beta %*% x))
        c3 = sum(lambda * (phi(alpha - beta %*% x)[2:(J + 1)] - phi(alpha -
                                                                      beta %*% x)[1:J]))
        gbln[j] = gbln[j] - b3 / c3
        gcln[j] = gcln[j] - a3 / c3
      }
    }
  }
  
  return(c(gbln, galn, gcln))
}




#the main  function to caculated the three estimathe and p-value
##the function use Rcpp function
cp.estimate = function(y, G, X) {
  n=table(y)
  J=length(n)
  N = sum(n)
  
  K = ifelse(is.null(X), 0, ncol(as.matrix(X)))
  #----------------- conbime likelihood-----------------
  #constraint
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
  mod = maxLik(
    likelihood,
    grad = likelihoodgrad,
    start = sta,
    constraints = list(ineqA = ui, ineqB = ci),
    y = y,
    X = as.matrix(cbind(G, X)),
    J = J,
    lambda1=n[1]/N
  )
  #the estimte coffiecent
  co = mod$estimate
  ccdata = data.frame(y = y, G = G, X = X)
  B = p.likelihood.fisher(co, data = df.to.list(ccdata))
  A = ginv(-mod$hessian / N)
  # the variance
  va =  A %*% B %*% A / N
  ##-------------------the weight----------------
  cp = co[3 * J + K - 1]
  cm = co[J]
  cpro.var = va[3 * J + K - 1, 3 * J + K - 1]
  cmod.var = va[J, J]
  rho = va[J, 3 * J + K - 1]/sqrt(cpro.var*cmod.var)
  
  #--------------------------------pro and mod---------------------
  pro.mod = lrm(y ~ ., data = ccdata)
  p = pro.mod$coefficients[J]
  pro.var=pro.mod$var[J,J]
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

  m = mod.model$estimate[J]
  # mod.var=ginv(-mod.model$hessian)[J,J]
  # mod.var=ifelse(is.na(mod.var),cmod.var,mod.var)
  mod.var=cmod.var
  mod.statistic = m / sqrt(mod.var)
  ##p-value
  pro.p = anova(pro.mod)[1, 3]
  mod.p = 2 * (1 - pnorm(abs(mod.statistic)))
  
  covar = rho*sqrt(pro.var*mod.var)
  
  
  ##the weight test statisitc
  cov = matrix(c(pro.var, covar, covar, mod.var), ncol = 2)
  a = pro.var + (p - m) ^ 2
  b = mod.var
  c = covar
  w= optimize(function(w,a,b,c){w^2*a+2*w*(1-w)*c+(1-w)^2*b},c(0,1),a=a,b=b,c=c)$minimum
  new=w*p+(1-w)*m
  #the wald test statistic
  cob = c(w, 1 - w)
  new.var = t(cob) %*% cov %*% cob
  new.statistic = new / sqrt(new.var)
  new.p = 2 * (1 - pnorm(abs(new.statistic)))
  
  
  out=c(p, m, new,pro.p,mod.p,new.p,rho)
  # names(out)=c("pro","mod","weight","cpT.p","cmT.p","wT.p")
  return(out)
  
}
