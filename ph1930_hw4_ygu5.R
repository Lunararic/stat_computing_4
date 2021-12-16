author="yue.gu@uth.tmc.edu"

## 1 EM
em_mix = function(y, maxit=1000, tol=1e-4){
  counter = 0
  myP = c(1/3, 1/3, 1/3)
  #cat("iter = ", counter, "P = ", myP, "\n")
  myMu1 = rnorm(3, mean(y), 10)
  nsub = length(y)
  
  # E-step
  while (counter <= maxit) {
    counter = counter + 1
    # given mixture distribution of 3 normal distributions with equal variance 1 and different means
    # update p
    denominator = (myP[1]*dnorm(Y, mean=myMu1[1], sd=sqrt(1)) 
                   + myP[2]*dnorm(Y, mean=myMu1[2], sd=sqrt(1)) 
                   + (1 - myP[1] - myP[2])*dnorm(Y, mean=myMu1[3], sd=sqrt(1)))
    
    myP_new = c(sum(myP[1]*dnorm(Y, mean=myMu1[1], sd=sqrt(1))/denominator)/nsub, 
                sum(myP[2]*dnorm(Y, mean=myMu1[2], sd=sqrt(1))/denominator)/nsub,
                sum((1 - myP[1] - myP[2])*dnorm(Y, mean=myMu1[3], sd=sqrt(1))/denominator)/nsub)
    
    
    #cat("iter = ", counter, "P = ", myP.new, "\n")
    
    # update mu
    for (i in seq(myMu1)) {
      nominator = sum(((myP[i] * dnorm(y, mean=myMu1[i], sd=1))/denominator) * y) 
      mu_new[i] = nominator/sum((myP[i] * dnorm(y, mean=myMu1[i], sd=1))/denominator)
    }
    
    

    # stop iteration when exceed tolerance
    diff = max(c(base::norm(matrix(myP.new - myP), type="F"), base::norm(matrix(na.omit(mu_new - myMu1)), "F")))
    if (diff < tol){
      return(list(P = myP))
    } else {
      myP = myP.new
      myMu1 = mu_new
    }
  }
    
    # M-step
  return(list(P=myP, mu=myMu1))
  
}


# test
# simulate data
# set.seed(10)
# myP <- c(1/3, 1/6, 1/2)
# myMu <- c(8, 12, 25)
# mySigmaSq <- 1
# nsub <- 2000
# Y <- NULL
# for (i in 1:nsub) {
#  token <- sample(c(1, 2, 3), size=1, prob=myP)
#  if (token == 1) {
#    Y <- c(Y, rnorm(1, mean = myMu[1], sd = 1))
#  } else if(token == 2) {
#    Y <- c(Y, rnorm(1, mean = myMu[2], sd = 1))
#  } else {
#    Y <- c(Y, rnorm(1, mean = myMu[3], sd = 1))
#  }
# }
# 
# plot(density(Y))
# y = Y
# em_mix(Y)


# lecture codes
#myEM <- function(Y, tol=1e-6, maxit=1000){
#  counter <- 0
#  myP <- 0.5
#  cat("iter = ", counter, "P = ", myP, "\n")
#  while (counter <= maxit) {
#    counter <- counter + 1
#    denominator <- (myP*dnorm(Y, mean=myMu1, sd=sqrt(mySigmaSq))
#                    + (1-myP)*dnorm(Y, mean=myMu2, sd=sqrt(mySigmaSq)))
#    myP.new <- sum(myP*dnorm(Y, mean=myMu1, sd=sqrt(mySigmaSq))/denominator)/nsub
#    cat("iter = ", counter, "P = ", myP.new, "\n")
#    
#    if (abs(myP.new-myP) < tol){  # calculate the l1 norm
#      cat("\nSuccessfully Converged\n")
#      cat("iter = ", counter, "P = ", myP.new, "\n")
#      return(list(P = myP))
#    } else {
#      myP <- myP.new
#    }
#  }
#  print("Convergence Failed")
#  return(list(P = myP))
#}  


## 2 Optimization
admm_elastic = function(X, y, lambda, maxit=1000, tol=1e-4){
  
  alpha = 0.95
  # set soft threshold
  softThresh = function(x, lambda) {
    sign(x) * pmax(0, abs(x) - lambda)
  }  
  
  # pre-computation
  XX = t(X) %*% X
  Xy = t(X) %*% y
  # number of parameters
  p = ncol(X)
  
  # Lagrange Multiplier
  lambda = rep(0, p)
  # initiate rho as penalty coefficient and set maximum value
  maxRho = 5
  rho = 4
  
  z0 = z = beta0 = beta = rep(0, p)
  # calculate inverse matrix
  Sinv = solve(XX + rho * diag(rep(1, p)) + lambda * (1 - alpha) * diag(rep(1, p)))
  
  
  for (it in 1:maxit) {
    ## update beta
    ## beta <- solve(XX + rho*diag(rep(1, p)) ) %*% (Xy + rho * z - lambda)
    beta = Sinv %*% (Xy + rho * z - lambda)
    
    ## update z
    z = softThresh(beta + (lambda/rho), (lambda*alpha)/rho)
    
    ## update lambda
    lambda = lambda + rho* (beta - z ) 
    ## increase rho
    ## rho <- min(maxRho, rho*1.1)
    change = max(c(base::norm(beta - beta0, "F"),
                   base::norm(z - z0, "F")))
    if (change < tol || it > maxit) {
      break
    }
    
    # save last iteration results
    beta0 = beta
    z0 = z
  }
  return(beta0)
}

# test
# simulate data
# alpha <- 0.95
# set.seed(10)
# n <- 100
# p <- 500
# X <- matrix(rnorm(n*p), n, p)
# b <- rep(0, 500)
# b[301:305] <- c(5:1)*5
# y <-  X%*%b + rnorm(n)
# ynew <-  X%*%b + rnorm(n)
# admm_elastic(X, y, n)

# lecture codes
#softThresh <- function(x, lambda) {
#  sign(x)*pmax(0, abs(x) - lambda)
#}  

#admmLasso <- function(X, y, tau, maxit = 1000, tol=1e-4) {
#  XX <- t(X) %*% X
#  Xy <- t(X) %*% y
#  p <- ncol(X)
#  lambda <- rep(0, p)
#  maxRho <- 5
#  rho <- 4
#  z0 <- z <- beta0 <- beta <- rep(0, p)
#  Sinv <- solve(XX + rho*diag(rep(1, p)) )
#  for (it in 1:maxit) {
#    ## update beta
#    ## beta <- solve(XX + rho*diag(rep(1, p)) ) %*% (Xy + rho * z - lambda)
#    beta <- Sinv %*% (Xy + rho * z - lambda)
#    ## update z
#    z <- softThresh(beta + lambda/rho, tau/rho)
#    ## update lambda
#    lambda <- lambda + rho* (beta - z ) 
#    ## increase rho
#    ## rho <- min(maxRho, rho*1.1)
#    change <- max(  c( base::norm(beta - beta0, "F"),
#                       base::norm(z - z0, "F") ) )
#    if (change < tol || it > maxit) {
#      break
#    }
#    beta0 <-  beta
#    z0 <-  z
#  }
#  return(z)
#}       


## 3 Improve Speed
# calling C++ to improve running speed
install.packages("Rcpp", repos="https://cloud.r-project.org")
install.packages("RcppArmadillo", repos="https://cloud.r-project.org")

library(Rcpp)
library(RcppArmadillo)

admm_elastic_fast = function(X, y, lambda, maxit=1000, tol=1e-4){
  sourceCpp("admm_lasso_elastic_net.cpp")
  return(admm_lasso_elastic_net(X, y, lambda, maxit=1000, tol=1e-4))
}

# test
# admm_elastic_fast(X, y, n)
