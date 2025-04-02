# #example 1 AJEL
library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS########################

n = 50
# rho = c(-0.9, -0.5, 0, 0.5, 0.75)
rho = c(-0.6, -0.5, 0, 0.5, 0.6)
iter = 100

########Define Matrix of Results###################

coverage.jel = vector()
length.jel = vector()
upper.jel = vector()
lower.jel = vector()
results = matrix(NA, nrow = 5, ncol = 4)

########## FUNCTIONS#########################

get.NWK = function(x,u,small.u){
  bw = density(U, kernel = c("epanechnikov"))$bw
  KX = vector()
  K = vector()
  for(j in 1:length(u)){
    d = (u[j]-small.u)/bw
    KX[j] = exp(x[j])*(bw^(-1)*max(0,0.75*(1-d^2)))
    K[j] = (bw^(-1)*max(0,0.75*(1-d^2)))
  }
  result = sum(KX)/sum(K)
  return(result)
}

findci <- function(x.vector, AJEL = FALSE) {
  cifunc <- function(ci.val, x.vector, AJEL = AJEL) {
    if (AJEL == TRUE) {
      x.vector = x.vector - ci.val
      x.vector[length(x.vector)+1] = -0.5*log(length(x.vector))*mean(x.vector)
      el.test(x.vector,0)
    } else {el.test(x.vector - ci.val,0)}
  }
  
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector, AJEL = AJEL)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector, AJEL = AJEL)$Up
  return(list(Up=U, Low = L))
}

############ITERATIONS#######################

for (ii in 1:5){
  for (jj in 1:iter){
    
    ####### DATA GENERATION#####################
    
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(-1,1), Sigma = sigma ) 
    colnames(data) <- c("X","Y")
    
    U = runif(n,0,1)
    
    ########DISTORTING FUNCTION#################
    
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    
    ######OBSERVED DATA WITH MULTIPLICATIVE DISTORTION###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    #######CALCULATING ESTIMATORS##############
    
    #RESIDUAL BASED ESTIMATOR###
    
    e.Y.est = vector()
    for (i in 1:n){
      e.Y.est[i] = xy.obs[,2][i] - log(get.NWK(xy.obs[,2],U,U[i])) + log(mean(exp(xy.obs[,2])))
    }
    
    e.X.est = vector()
    for (i in 1:n){
      e.X.est[i] = xy.obs[,1][i] - log(get.NWK(xy.obs[,1],U,U[i])) + log(mean(exp(xy.obs[,1])))
    }
    
    cov.e.est = mean(e.X.est*e.Y.est) - mean(e.X.est)*mean(e.Y.est)
    sig.e.X = mean(e.X.est*e.X.est) - mean(e.X.est)*mean(e.X.est)
    sig.e.Y = mean(e.Y.est*e.Y.est) - mean(e.Y.est)*mean(e.Y.est)
    
    rho.e.est = (cov.e.est)/(sqrt(sig.e.X*sig.e.Y))
    
    ############JACKKNIFING#######################
    
    rho.j=vector()
    
    for (j in 1:n){
      
      xy.jack=xy.obs[-j,]
      U.jack = U[-j]
      
      Eix.j = vector()
      Eiy.j = vector()
      
      for (k in 1:(n-1)){
        
        Eix.j[k] = xy.jack[,1][k]-log(get.NWK(xy.jack[,1],U.jack,U.jack[k]))+log(mean(exp(xy.jack[,1])))
        Eiy.j[k] =xy.jack[,2][k]- log(get.NWK(xy.jack[,2],U.jack,U.jack[k]))+log(mean(exp(xy.jack[,2])))
      }
      
      cov.e.est.jack = mean(Eix.j*Eiy.j) - mean(Eix.j)*mean(Eiy.j)
      sig.e.X.jack = mean(Eix.j*Eix.j) - mean(Eix.j)*mean(Eix.j)
      sig.e.Y.jack = mean(Eiy.j*Eiy.j) - mean(Eiy.j)*mean(Eiy.j)
      
      rho.e.est.jack = (cov.e.est.jack)/(sqrt(sig.e.X.jack*sig.e.Y.jack))
      
      rho.j[j]=n*rho.e.est-(n-1)*rho.e.est.jack
    }
    
    #########JACKKNIFE EMPIRICAL LIKELIHOOD#################
    
    wni = rho.j - rho[ii]
    wni.ci = rho.j - mean(rho.j)
    
    il.jeltest = el.test(wni,0)
    il.j = il.jeltest$ '-2LLR'
    ci.jel = findci(rho.j)
    upper.jel[jj] = ci.jel$Up
    lower.jel[jj] = ci.jel$Low
    length.jel[jj] = upper.jel[jj] - lower.jel[jj]
    coverage.jel[jj] = il.j < 1.96^2
  }
  
  results[ii,1] = mean(lower.jel)
  results[ii,2] = mean(upper.jel)
  results[ii,3] = mean(length.jel)
  results[ii,4] = sum(coverage.jel)/iter
}

results





# 
# [,1]       [,2]      [,3] [,4]
# [1,] -0.8795264 -0.6448200 0.2347064 0.36
# [2,] -0.6615046 -0.2116396 0.4498650 0.95
# [3,] -0.3339047  0.2553070 0.5892116 0.96
# [4,]  0.0726691  0.5944256 0.5217565 0.79
# [5,]  0.3292181  0.7406511 0.4114329 0.47




# #example 1 EL

library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS########################

n = 50
# rho = c(-0.9, -0.5, 0, 0.5, 0.75)
rho = c(-0.6, -0.5, 0, 0.5, 0.6)
iter = 100

########Define Matrix of Results###################

coverage.jel = vector()
length.jel = vector()
upper.jel = vector()
lower.jel = vector()
results = matrix(NA, nrow = 5, ncol = 4)

########## FUNCTIONS#########################

get.NWK = function(x,u,small.u){
  bw = density(U, kernel = c("epanechnikov"))$bw
  KX = vector()
  K = vector()
  for(j in 1:length(u)){
    d = (u[j]-small.u)/bw
    KX[j] = exp(x[j])*(bw^(-1)*max(0,0.75*(1-d^2)))
    K[j] = (bw^(-1)*max(0,0.75*(1-d^2)))
  }
  result = sum(KX)/sum(K)
  return(result)
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    el.test(x.vector - ci.val,0)
  }
  
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up=U, Low = L))
}

############ITERATIONS#######################

for (ii in 1:5){
  for (jj in 1:iter){
    
    ####### DATA GENERATION#####################
    
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(-1,1), Sigma = sigma ) 
    colnames(data) <- c("X","Y")
    
    U = runif(n,0,1)
    
    ########DISTORTING FUNCTION#################
    
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    
    ######OBSERVED DATA WITH MULTIPLICATIVE DISTORTION###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    #######CALCULATING ESTIMATORS##############
    
    #RESIDUAL BASED ESTIMATOR###
    
    e.Y.est = vector()
    for (i in 1:n){
      e.Y.est[i] = xy.obs[,2][i] - log(get.NWK(xy.obs[,2],U,U[i])) + log(mean(exp(xy.obs[,2])))
    }
    
    e.X.est = vector()
    for (i in 1:n){
      e.X.est[i] = xy.obs[,1][i] - log(get.NWK(xy.obs[,1],U,U[i])) + log(mean(exp(xy.obs[,1])))
    }
    
    cov.e.est = mean(e.X.est*e.Y.est) - mean(e.X.est)*mean(e.Y.est)
    sig.e.X = mean(e.X.est*e.X.est) - mean(e.X.est)*mean(e.X.est)
    sig.e.Y = mean(e.Y.est*e.Y.est) - mean(e.Y.est)*mean(e.Y.est)
    
    rho.e.est = (cov.e.est)/(sqrt(sig.e.X*sig.e.Y))
    
    ############JACKKNIFING#######################
    
    rho.j=vector()
    
    for (j in 1:n){
      
      xy.jack=xy.obs[-j,]
      U.jack = U[-j]
      
      Eix.j = vector()
      Eiy.j = vector()
      
      for (k in 1:(n-1)){
        
        Eix.j[k] = xy.jack[,1][k]-log(get.NWK(xy.jack[,1],U.jack,U.jack[k]))+log(mean(exp(xy.jack[,1])))
        Eiy.j[k] =xy.jack[,2][k]- log(get.NWK(xy.jack[,2],U.jack,U.jack[k]))+log(mean(exp(xy.jack[,2])))
      }
      
      cov.e.est.jack = mean(Eix.j*Eiy.j) - mean(Eix.j)*mean(Eiy.j)
      sig.e.X.jack = mean(Eix.j*Eix.j) - mean(Eix.j)*mean(Eix.j)
      sig.e.Y.jack = mean(Eiy.j*Eiy.j) - mean(Eiy.j)*mean(Eiy.j)
      
      rho.e.est.jack = (cov.e.est.jack)/(sqrt(sig.e.X.jack*sig.e.Y.jack))
      
      rho.j[j]=n*rho.e.est-(n-1)*rho.e.est.jack
    }
    
    #########JACKKNIFE EMPIRICAL LIKELIHOOD#################
    
    wni = rho.j - rho[ii]
    wni.ci = rho.j - mean(rho.j)
    
    il.jeltest = el.test(wni,0)
    il.j = il.jeltest$ '-2LLR'
    ci.jel = findci(rho.j)
    upper.jel[jj] = ci.jel$Up
    lower.jel[jj] = ci.jel$Low
    length.jel[jj] = upper.jel[jj] - lower.jel[jj]
    coverage.jel[jj] = il.j < 1.96^2
  }
  
  results[ii,1] = mean(lower.jel)
  results[ii,2] = mean(upper.jel)
  results[ii,3] = mean(length.jel)
  results[ii,4] = sum(coverage.jel)/iter
}

results


# 
# [,1]       [,2]      [,3] [,4]
# [1,] -0.86852548 -0.6338624 0.2346630 0.32
# [2,] -0.66801863 -0.2064292 0.4615895 0.91
# [3,] -0.33661268  0.2387321 0.5753447 0.95
# [4,]  0.08099838  0.6046639 0.5236655 0.78
# [5,]  0.34022734  0.7613105 0.4210831 0.46






# #example 1 JEL
library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS########################

n = 50
# rho = c(-0.9, -0.5, 0, 0.5, 0.75)
rho = c(-0.6, -0.5, 0, 0.5, 0.6)
iter = 100

########Define Matrix of Results###################

coverage.jel = vector()
length.jel = vector()
upper.jel = vector()
lower.jel = vector()
results = matrix(NA, nrow = 5, ncol = 4)

########## FUNCTIONS#########################

get.NWK = function(x,u,small.u){
  bw = density(U, kernel = c("epanechnikov"))$bw
  KX = vector()
  K = vector()
  for(j in 1:length(u)){
    d = (u[j]-small.u)/bw
    KX[j] = exp(x[j])*(bw^(-1)*max(0,0.75*(1-d^2)))
    K[j] = (bw^(-1)*max(0,0.75*(1-d^2)))
  }
  result = sum(KX)/sum(K)
  return(result)
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    el.test(x.vector - ci.val,0)
  }
  
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up=U, Low = L))
}

############ITERATIONS#######################

for (ii in 1:5){
  for (jj in 1:iter){
    
    ####### DATA GENERATION#####################
    
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(-1,1), Sigma = sigma ) 
    colnames(data) <- c("X","Y")
    
    U = runif(n,0,1)
    
    ########DISTORTING FUNCTION#################
    
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    
    ######OBSERVED DATA WITH MULTIPLICATIVE DISTORTION###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    #######CALCULATING ESTIMATORS##############
    
    #RESIDUAL BASED ESTIMATOR###
    
    e.Y.est = vector()
    for (i in 1:n){
      e.Y.est[i] = xy.obs[,2][i] - log(get.NWK(xy.obs[,2],U,U[i])) + log(mean(exp(xy.obs[,2])))
    }
    
    e.X.est = vector()
    for (i in 1:n){
      e.X.est[i] = xy.obs[,1][i] - log(get.NWK(xy.obs[,1],U,U[i])) + log(mean(exp(xy.obs[,1])))
    }
    
    cov.e.est = mean(e.X.est*e.Y.est) - mean(e.X.est)*mean(e.Y.est)
    sig.e.X = mean(e.X.est*e.X.est) - mean(e.X.est)*mean(e.X.est)
    sig.e.Y = mean(e.Y.est*e.Y.est) - mean(e.Y.est)*mean(e.Y.est)
    
    rho.e.est = (cov.e.est)/(sqrt(sig.e.X*sig.e.Y))
    
    ############JACKKNIFING#######################
    
    rho.j=vector()
    
    for (j in 1:n){
      
      xy.jack=xy.obs[-j,]
      U.jack = U[-j]
      
      Eix.j = vector()
      Eiy.j = vector()
      
      for (k in 1:(n-1)){
        
        Eix.j[k] = xy.jack[,1][k]-log(get.NWK(xy.jack[,1],U.jack,U.jack[k]))+log(mean(exp(xy.jack[,1])))
        Eiy.j[k] =xy.jack[,2][k]- log(get.NWK(xy.jack[,2],U.jack,U.jack[k]))+log(mean(exp(xy.jack[,2])))
      }
      
      cov.e.est.jack = mean(Eix.j*Eiy.j) - mean(Eix.j)*mean(Eiy.j)
      sig.e.X.jack = mean(Eix.j*Eix.j) - mean(Eix.j)*mean(Eix.j)
      sig.e.Y.jack = mean(Eiy.j*Eiy.j) - mean(Eiy.j)*mean(Eiy.j)
      
      rho.e.est.jack = (cov.e.est.jack)/(sqrt(sig.e.X.jack*sig.e.Y.jack))
      
      rho.j[j]=n*rho.e.est-(n-1)*rho.e.est.jack
    }
    
    #########JACKKNIFE EMPIRICAL LIKELIHOOD (JEL) #################
    
    wni = rho.j - rho[ii]
    wni.ci = rho.j - mean(rho.j)
    
    il.jeltest = el.test(wni,0)
    il.j = il.jeltest$ '-2LLR'
    ci.jel = findci(rho.j)
    upper.jel[jj] = ci.jel$Up
    lower.jel[jj] = ci.jel$Low
    length.jel[jj] = upper.jel[jj] - lower.jel[jj]
    coverage.jel[jj] = il.j < 1.96^2
  }
  
  results[ii,1] = mean(lower.jel)
  results[ii,2] = mean(upper.jel)
  results[ii,3] = mean(length.jel)
  results[ii,4] = sum(coverage.jel)/iter
}

results

# 
# [,1]       [,2]      [,3] [,4]
# [1,] -0.87648623 -0.6439523 0.2325339 0.24
# [2,] -0.66943360 -0.2372104 0.4322232 0.93
# [3,] -0.37192650  0.2195131 0.5914396 0.82
# [4,]  0.07136188  0.6007213 0.5293594 0.80
# [5,]  0.35388200  0.7679637 0.4140817 0.61






# #example 1 MJEL


library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS########################

n = 50
# rho = c(-0.9, -0.5, 0, 0.5, 0.75)
rho = c(-0.6, -0.5, 0, 0.5, 0.6)
iter = 100

########Define Matrix of Results###################

coverage.jel = vector()
length.jel = vector()
upper.jel = vector()
lower.jel = vector()
results = matrix(NA, nrow = 5, ncol = 4)

########## FUNCTIONS#########################

get.NWK = function(x,u,small.u){
  bw = density(U, kernel = c("epanechnikov"))$bw
  KX = vector()
  K = vector()
  for(j in 1:length(u)){
    d = (u[j]-small.u)/bw
    KX[j] = exp(x[j])*(bw^(-1)*max(0,0.75*(1-d^2)))
    K[j] = (bw^(-1)*max(0,0.75*(1-d^2)))
  }
  result = sum(KX)/sum(K)
  return(result)
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    n = length(x.vector)
    adjusted_mean = mean(x.vector) - 0.5 * log(n) * mean(x.vector)
    x.vector = x.vector - ci.val
    x.vector = c(x.vector, adjusted_mean)
    el.test(x.vector, 0)
  }
  
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up=U, Low = L))
}

############ITERATIONS#######################

for (ii in 1:5){
  for (jj in 1:iter){
    
    ####### DATA GENERATION#####################
    
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(-1,1), Sigma = sigma ) 
    colnames(data) <- c("X","Y")
    
    U = runif(n,0,1)
    
    ########DISTORTING FUNCTION#################
    
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    
    ######OBSERVED DATA WITH MULTIPLICATIVE DISTORTION###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    #######CALCULATING ESTIMATORS##############
    
    #RESIDUAL BASED ESTIMATOR###
    
    e.Y.est = vector()
    for (i in 1:n){
      e.Y.est[i] = xy.obs[,2][i] - log(get.NWK(xy.obs[,2],U,U[i])) + log(mean(exp(xy.obs[,2])))
    }
    
    e.X.est = vector()
    for (i in 1:n){
      e.X.est[i] = xy.obs[,1][i] - log(get.NWK(xy.obs[,1],U,U[i])) + log(mean(exp(xy.obs[,1])))
    }
    
    cov.e.est = mean(e.X.est*e.Y.est) - mean(e.X.est)*mean(e.Y.est)
    sig.e.X = mean(e.X.est*e.X.est) - mean(e.X.est)*mean(e.X.est)
    sig.e.Y = mean(e.Y.est*e.Y.est) - mean(e.Y.est)*mean(e.Y.est)
    
    rho.e.est = (cov.e.est)/(sqrt(sig.e.X*sig.e.Y))
    
    ############JACKKNIFING#######################
    
    rho.j=vector()
    
    for (j in 1:n){
      
      xy.jack=xy.obs[-j,]
      U.jack = U[-j]
      
      Eix.j = vector()
      Eiy.j = vector()
      
      for (k in 1:(n-1)){
        
        Eix.j[k] = xy.jack[,1][k]-log(get.NWK(xy.jack[,1],U.jack,U.jack[k]))+log(mean(exp(xy.jack[,1])))
        Eiy.j[k] =xy.jack[,2][k]- log(get.NWK(xy.jack[,2],U.jack,U.jack[k]))+log(mean(exp(xy.jack[,2])))
      }
      
      cov.e.est.jack = mean(Eix.j*Eiy.j) - mean(Eix.j)*mean(Eiy.j)
      sig.e.X.jack = mean(Eix.j*Eix.j) - mean(Eix.j)*mean(Eix.j)
      sig.e.Y.jack = mean(Eiy.j*Eiy.j) - mean(Eiy.j)*mean(Eiy.j)
      
      rho.e.est.jack = (cov.e.est.jack)/(sqrt(sig.e.X.jack*sig.e.Y.jack))
      
      rho.j[j]=n*rho.e.est-(n-1)*rho.e.est.jack
    }
    
    #########JACKKNIFE MODIFIED EMPIRICAL LIKELIHOOD (MJEL) #################
    
    wni = rho.j - rho[ii]
    wni.ci = rho.j - mean(rho.j)
    
    il.jeltest = el.test(wni,0)
    il.j = il.jeltest$ '-2LLR'
    ci.jel = findci(rho.j)
    upper.jel[jj] = ci.jel$Up
    lower.jel[jj] = ci.jel$Low
    length.jel[jj] = upper.jel[jj] - lower.jel[jj]
    coverage.jel[jj] = il.j < 1.96^2
  }
  
  results[ii,1] = mean(lower.jel)
  results[ii,2] = mean(upper.jel)
  results[ii,3] = mean(length.jel)
  results[ii,4] = sum(coverage.jel)/iter
}

results


# [,1]       [,2]      [,3] [,4]
# [1,] -0.85899341 -0.6097990 0.2491945 0.28
# [2,] -0.66532107 -0.2196707 0.4456504 0.94
# [3,] -0.31635689  0.2705197 0.5868766 0.95
# [4,]  0.07333963  0.5908701 0.5175305 0.78
# [5,]  0.28767708  0.7289504 0.4412733 0.46





# #example 1 TJEL

library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS########################

n = 50
# rho = c(-0.9, -0.5, 0, 0.5, 0.75)
rho = c(-0.6, -0.5, 0, 0.5, 0.6)
iter = 100

########Define Matrix of Results###################

coverage.jel = vector()
length.jel = vector()
upper.jel = vector()
lower.jel = vector()
results = matrix(NA, nrow = 5, ncol = 4)

########## FUNCTIONS#########################

get.NWK = function(x,u,small.u){
  bw = density(U, kernel = c("epanechnikov"))$bw
  KX = vector()
  K = vector()
  for(j in 1:length(u)){
    d = (u[j]-small.u)/bw
    KX[j] = exp(x[j])*(bw^(-1)*max(0,0.75*(1-d^2)))
    K[j] = (bw^(-1)*max(0,0.75*(1-d^2)))
  }
  result = sum(KX)/sum(K)
  return(result)
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    n = length(x.vector)
    tau = sqrt(n) * mean(x.vector)
    x.vector = x.vector - ci.val
    x.vector = c(x.vector, tau)
    el.test(x.vector, 0)
  }
  
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up=U, Low = L))
}

############ITERATIONS#######################

for (ii in 1:5){
  for (jj in 1:iter){
    
    ####### DATA GENERATION#####################
    
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(-1,1), Sigma = sigma ) 
    colnames(data) <- c("X","Y")
    
    U = runif(n,0,1)
    
    ########DISTORTING FUNCTION#################
    
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    
    ######OBSERVED DATA WITH MULTIPLICATIVE DISTORTION###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    #######CALCULATING ESTIMATORS##############
    
    #RESIDUAL BASED ESTIMATOR###
    
    e.Y.est = vector()
    for (i in 1:n){
      e.Y.est[i] = xy.obs[,2][i] - log(get.NWK(xy.obs[,2],U,U[i])) + log(mean(exp(xy.obs[,2])))
    }
    
    e.X.est = vector()
    for (i in 1:n){
      e.X.est[i] = xy.obs[,1][i] - log(get.NWK(xy.obs[,1],U,U[i])) + log(mean(exp(xy.obs[,1])))
    }
    
    cov.e.est = mean(e.X.est*e.Y.est) - mean(e.X.est)*mean(e.Y.est)
    sig.e.X = mean(e.X.est*e.X.est) - mean(e.X.est)*mean(e.X.est)
    sig.e.Y = mean(e.Y.est*e.Y.est) - mean(e.Y.est)*mean(e.Y.est)
    
    rho.e.est = (cov.e.est)/(sqrt(sig.e.X*sig.e.Y))
    
    ############JACKKNIFING#######################
    
    rho.j=vector()
    
    for (j in 1:n){
      
      xy.jack=xy.obs[-j,]
      U.jack = U[-j]
      
      Eix.j = vector()
      Eiy.j = vector()
      
      for (k in 1:(n-1)){
        
        Eix.j[k] = xy.jack[,1][k]-log(get.NWK(xy.jack[,1],U.jack,U.jack[k]))+log(mean(exp(xy.jack[,1])))
        Eiy.j[k] =xy.jack[,2][k]- log(get.NWK(xy.jack[,2],U.jack,U.jack[k]))+log(mean(exp(xy.jack[,2])))
      }
      
      cov.e.est.jack = mean(Eix.j*Eiy.j) - mean(Eix.j)*mean(Eiy.j)
      sig.e.X.jack = mean(Eix.j*Eix.j) - mean(Eix.j)*mean(Eix.j)
      sig.e.Y.jack = mean(Eiy.j*Eiy.j) - mean(Eiy.j)*mean(Eiy.j)
      
      rho.e.est.jack = (cov.e.est.jack)/(sqrt(sig.e.X.jack*sig.e.Y.jack))
      
      rho.j[j]=n*rho.e.est-(n-1)*rho.e.est.jack
    }
    
    #########JACKKNIFE TRANSFORMED EMPIRICAL LIKELIHOOD (TJEL) #################
    
    wni = rho.j - rho[ii]
    wni.ci = rho.j - mean(rho.j)
    
    il.jeltest = el.test(wni,0)
    il.j = il.jeltest$ '-2LLR'
    ci.jel = findci(rho.j)
    upper.jel[jj] = ci.jel$Up
    lower.jel[jj] = ci.jel$Low
    length.jel[jj] = upper.jel[jj] - lower.jel[jj]
    coverage.jel[jj] = il.j < 1.96^2
  }
  
  results[ii,1] = mean(lower.jel)
  results[ii,2] = mean(upper.jel)
  results[ii,3] = mean(length.jel)
  results[ii,4] = sum(coverage.jel)/iter
}

results
# 
# 
# [,1]       [,2]      [,3] [,4]
# [1,] -1.2569518 -0.6948707 0.5620812 0.27
# [2,] -0.7934466 -0.2580921 0.5353545 0.93
# [3,] -0.3457428  0.2526174 0.5983602 0.96
# [4,]  0.1316587  0.7162662 0.5846075 0.85
# [5,]  0.3977943  0.9552846 0.5574903 0.49

