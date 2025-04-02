# #example 1 AJEL with Beta-distributed data
# install.packages('emplik')
# install.packages('kedd')

library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS########################

n = 75
rho = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.6)
iter = 100

########Define Matrix of Results###################

coverage.jel = vector()
length.jel = vector()
upper.jel = vector()
lower.jel = vector()
results = matrix(NA, nrow = length(rho), ncol = 4)

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

for (ii in 1:length(rho)){
  
  for (jj in 1:iter){
    
    ####### DATA GENERATION using BETA DISTRIBUTION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    norm_data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
    
    U1 <- pnorm(norm_data[, 1])
    U2 <- pnorm(norm_data[, 2])
    
    X <- qbeta(U1, shape1 = 2, shape2 = 5)
    Y <- qbeta(U2, shape1 = 2, shape2 = 5)
    
    data <- cbind(X, Y)
    colnames(data) <- c("X", "Y")
    
    U = runif(n,0,1)
    
    ########DISTORTING FUNCTION#################
    
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    
    ######OBSERVED DATA WITH MULTIPLICATIVE DISTORTION###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    #######CALCULATING ESTIMATORS##############
    
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




# #example 1 EL with Beta-distributed data
# install.packages('emplik')
# install.packages('kedd')

library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS########################

n = 75
rho = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.6)
iter = 100

coverage.el = vector()
length.el = vector()
upper.el = vector()
lower.el = vector()
results = matrix(NA, nrow = length(rho), ncol = 4)

########## FUNCTIONS #########################

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
    el.test(x.vector - ci.val, 0)
  }
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up=U, Low = L))
}

############ ITERATIONS #######################

for (ii in 1:length(rho)){
  
  for (jj in 1:iter){
    
    ####### DATA GENERATION using BETA DISTRIBUTION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    norm_data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
    
    U1 <- pnorm(norm_data[, 1])
    U2 <- pnorm(norm_data[, 2])
    
    X <- qbeta(U1, shape1 = 2, shape2 = 5)
    Y <- qbeta(U2, shape1 = 2, shape2 = 5)
    
    data <- cbind(X, Y)
    colnames(data) <- c("X", "Y")
    
    U = runif(n,0,1)
    
    ########DISTORTING FUNCTION#################
    
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    
    ######OBSERVED DATA WITH MULTIPLICATIVE DISTORTION###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    #######CALCULATING RESIDUAL-BASED ESTIMATOR##############
    
    e.Y.est = vector()
    for (i in 1:n){
      e.Y.est[i] = xy.obs[,2][i] - log(get.NWK(xy.obs[,2],U,U[i])) + log(mean(exp(xy.obs[,2])))
    }
    
    e.X.est = vector()
    for (i in 1:n){
      e.X.est[i] = xy.obs[,1][i] - log(get.NWK(xy.obs[,1],U,U[i])) + log(mean(exp(xy.obs[,1])))
    }
    
    cov.e.est = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
    sig.e.X = mean(e.X.est^2) - mean(e.X.est)^2
    sig.e.Y = mean(e.Y.est^2) - mean(e.Y.est)^2
    
    rho.e.est = cov.e.est / sqrt(sig.e.X * sig.e.Y)
    
    ######## EMPIRICAL LIKELIHOOD ####################
    
    wni = (e.X.est - mean(e.X.est)) * (e.Y.est - mean(e.Y.est)) / sqrt(sig.e.X * sig.e.Y)
    wni = wni - rho[ii]  # Centering for EL
    
    il.eltest = el.test(wni, 0)
    il.el = il.eltest$'-2LLR'
    ci.el = findci(wni)
    
    upper.el[jj] = ci.el$Up
    lower.el[jj] = ci.el$Low
    length.el[jj] = upper.el[jj] - lower.el[jj]
    coverage.el[jj] = il.el < 1.96^2
  }
  
  results[ii,1] = mean(lower.el)
  results[ii,2] = mean(upper.el)
  results[ii,3] = mean(length.el)
  results[ii,4] = sum(coverage.el)/iter
}

results




# #example 1 JEL with Beta-distributed data
# install.packages('emplik')
# install.packages('kedd')

library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 75
rho = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.6)
iter = 100

coverage.jel = vector()
length.jel = vector()
upper.jel = vector()
lower.jel = vector()
results = matrix(NA, nrow = length(rho), ncol = 4)

########## FUNCTIONS #########################

get.NWK = function(x, u, small.u) {
  bw = density(U, kernel = "epanechnikov")$bw
  KX = vector()
  K = vector()
  for(j in 1:length(u)){
    d = (u[j] - small.u)/bw
    KX[j] = exp(x[j]) * (bw^(-1) * max(0, 0.75 * (1 - d^2)))
    K[j] = (bw^(-1) * max(0, 0.75 * (1 - d^2)))
  }
  result = sum(KX)/sum(K)
  return(result)
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    el.test(x.vector - ci.val, 0)
  }
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up=U, Low = L))
}

############ ITERATIONS #######################

for (ii in 1:length(rho)) {
  
  for (jj in 1:iter) {
    
    ####### DATA GENERATION using BETA DISTRIBUTION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    norm_data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
    
    U1 <- pnorm(norm_data[, 1])
    U2 <- pnorm(norm_data[, 2])
    
    X <- qbeta(U1, shape1 = 2, shape2 = 5)
    Y <- qbeta(U2, shape1 = 2, shape2 = 5)
    
    data <- cbind(X, Y)
    colnames(data) <- c("X", "Y")
    
    U = runif(n, 0, 1)
    
    ########DISTORTING FUNCTION #################
    
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    
    ###### OBSERVED DATA WITH MULTIPLICATIVE DISTORTION ###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    ####### CALCULATING ESTIMATORS ####################
    
    e.Y.est = vector()
    for (i in 1:n){
      e.Y.est[i] = xy.obs[,2][i] - log(get.NWK(xy.obs[,2], U, U[i])) + log(mean(exp(xy.obs[,2])))
    }
    
    e.X.est = vector()
    for (i in 1:n){
      e.X.est[i] = xy.obs[,1][i] - log(get.NWK(xy.obs[,1], U, U[i])) + log(mean(exp(xy.obs[,1])))
    }
    
    cov.e.est = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
    sig.e.X = mean(e.X.est^2) - mean(e.X.est)^2
    sig.e.Y = mean(e.Y.est^2) - mean(e.Y.est)^2
    
    rho.e.est = cov.e.est / sqrt(sig.e.X * sig.e.Y)
    
    ########## JACKKNIFING ###########################
    
    rho.j = vector()
    
    for (j in 1:n) {
      
      xy.jack = xy.obs[-j, ]
      U.jack = U[-j]
      
      Eix.j = vector()
      Eiy.j = vector()
      
      for (k in 1:(n - 1)) {
        Eix.j[k] = xy.jack[,1][k] - log(get.NWK(xy.jack[,1], U.jack, U.jack[k])) + log(mean(exp(xy.jack[,1])))
        Eiy.j[k] = xy.jack[,2][k] - log(get.NWK(xy.jack[,2], U.jack, U.jack[k])) + log(mean(exp(xy.jack[,2])))
      }
      
      cov.e.est.jack = mean(Eix.j * Eiy.j) - mean(Eix.j) * mean(Eiy.j)
      sig.e.X.jack = mean(Eix.j^2) - mean(Eix.j)^2
      sig.e.Y.jack = mean(Eiy.j^2) - mean(Eiy.j)^2
      
      rho.e.est.jack = cov.e.est.jack / sqrt(sig.e.X.jack * sig.e.Y.jack)
      
      rho.j[j] = n * rho.e.est - (n - 1) * rho.e.est.jack
    }
    
    ########## JACKKNIFE EMPIRICAL LIKELIHOOD ##################
    
    wni = rho.j - rho[ii]
    
    il.jeltest = el.test(wni, 0)
    il.j = il.jeltest$'-2LLR'
    
    ci.jel = findci(rho.j)
    upper.jel[jj] = ci.jel$Up
    lower.jel[jj] = ci.jel$Low
    length.jel[jj] = upper.jel[jj] - lower.jel[jj]
    coverage.jel[jj] = il.j < 1.96^2
  }
  
  results[ii, 1] = mean(lower.jel)
  results[ii, 2] = mean(upper.jel)
  results[ii, 3] = mean(length.jel)
  results[ii, 4] = sum(coverage.jel) / iter
}

results




# #example 1 MJEL with Beta-distributed data
# install.packages('emplik')
# install.packages('kedd')

library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 75
rho = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.6)
iter = 100

coverage.mjel = vector()
length.mjel = vector()
upper.mjel = vector()
lower.mjel = vector()
results = matrix(NA, nrow = length(rho), ncol = 4)

########## FUNCTIONS #########################

get.NWK = function(x, u, small.u) {
  bw = density(U, kernel = "epanechnikov")$bw
  KX = vector()
  K = vector()
  for(j in 1:length(u)){
    d = (u[j] - small.u)/bw
    KX[j] = exp(x[j]) * (bw^(-1) * max(0, 0.75 * (1 - d^2)))
    K[j] = (bw^(-1) * max(0, 0.75 * (1 - d^2)))
  }
  result = sum(KX)/sum(K)
  return(result)
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    el.test(x.vector - ci.val, 0)
  }
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up=U, Low = L))
}

############ ITERATIONS #######################

for (ii in 1:length(rho)) {
  
  for (jj in 1:iter) {
    
    ####### DATA GENERATION using BETA DISTRIBUTION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    norm_data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
    
    U1 <- pnorm(norm_data[, 1])
    U2 <- pnorm(norm_data[, 2])
    
    X <- qbeta(U1, shape1 = 2, shape2 = 5)
    Y <- qbeta(U2, shape1 = 2, shape2 = 5)
    
    data <- cbind(X, Y)
    colnames(data) <- c("X", "Y")
    
    U = runif(n, 0, 1)
    
    ########DISTORTING FUNCTION #################
    
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    
    ###### OBSERVED DATA WITH MULTIPLICATIVE DISTORTION ###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    ####### CALCULATING ESTIMATORS ####################
    
    e.Y.est = vector()
    for (i in 1:n){
      e.Y.est[i] = xy.obs[,2][i] - log(get.NWK(xy.obs[,2], U, U[i])) + log(mean(exp(xy.obs[,2])))
    }
    
    e.X.est = vector()
    for (i in 1:n){
      e.X.est[i] = xy.obs[,1][i] - log(get.NWK(xy.obs[,1], U, U[i])) + log(mean(exp(xy.obs[,1])))
    }
    
    cov.e.est = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
    sig.e.X = mean(e.X.est^2) - mean(e.X.est)^2
    sig.e.Y = mean(e.Y.est^2) - mean(e.Y.est)^2
    
    rho.e.est = cov.e.est / sqrt(sig.e.X * sig.e.Y)
    
    ########## JACKKNIFING ###########################
    
    rho.j = vector()
    
    for (j in 1:n) {
      
      xy.jack = xy.obs[-j, ]
      U.jack = U[-j]
      
      Eix.j = vector()
      Eiy.j = vector()
      
      for (k in 1:(n - 1)) {
        Eix.j[k] = xy.jack[,1][k] - log(get.NWK(xy.jack[,1], U.jack, U.jack[k])) + log(mean(exp(xy.jack[,1])))
        Eiy.j[k] = xy.jack[,2][k] - log(get.NWK(xy.jack[,2], U.jack, U.jack[k])) + log(mean(exp(xy.jack[,2])))
      }
      
      cov.e.est.jack = mean(Eix.j * Eiy.j) - mean(Eix.j) * mean(Eiy.j)
      sig.e.X.jack = mean(Eix.j^2) - mean(Eix.j)^2
      sig.e.Y.jack = mean(Eiy.j^2) - mean(Eiy.j)^2
      
      rho.e.est.jack = cov.e.est.jack / sqrt(sig.e.X.jack * sig.e.Y.jack)
      
      rho.j[j] = n * rho.e.est - (n - 1) * rho.e.est.jack
    }
    
    ########## MJEL (Mean-adjusted) EMPIRICAL LIKELIHOOD ##################
    
    rho.bar = mean(rho.j)
    wni = rho.j - rho.bar
    
    il.mjeltest = el.test(wni, 0)
    il.mjel = il.mjeltest$'-2LLR'
    
    ci.mjel = findci(rho.j)
    upper.mjel[jj] = ci.mjel$Up
    lower.mjel[jj] = ci.mjel$Low
    length.mjel[jj] = upper.mjel[jj] - lower.mjel[jj]
    coverage.mjel[jj] = il.mjel < 1.96^2
  }
  
  results[ii, 1] = mean(lower.mjel)
  results[ii, 2] = mean(upper.mjel)
  results[ii, 3] = mean(length.mjel)
  results[ii, 4] = sum(coverage.mjel) / iter
}

results




# TJEL Method with Beta-distributed data
# install.packages('emplik')
# install.packages('kedd')

library(MASS)
library(emplik)
library(kedd)

# Parameters
n = 75
rho = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.6)
iter = 100

# Results storage
coverage.tjel = vector()
length.tjel = vector()
upper.tjel = vector()
lower.tjel = vector()
results = matrix(NA, nrow = length(rho), ncol = 4)

# Kernel smoother function
get.NWK = function(x, u, small.u) {
  bw = density(U, kernel = "epanechnikov")$bw
  KX = K = vector()
  for (j in 1:length(u)) {
    d = (u[j] - small.u) / bw
    KX[j] = exp(x[j]) * (bw^(-1) * max(0, 0.75 * (1 - d^2)))
    K[j] = (bw^(-1) * max(0, 0.75 * (1 - d^2)))
  }
  return(sum(KX) / sum(K))
}

# CI finder for transformed jackknife values
findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    el.test(x.vector - ci.val, 0)
  }
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up = U, Low = L))
}

# Simulation loop
for (ii in 1:length(rho)) {
  for (jj in 1:iter) {
    
    # Generate correlated beta data
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), 2, 2)
    norm_data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
    U1 <- pnorm(norm_data[, 1])
    U2 <- pnorm(norm_data[, 2])
    X <- qbeta(U1, 2, 5)
    Y <- qbeta(U2, 2, 5)
    data <- cbind(X, Y)
    colnames(data) <- c("X", "Y")
    
    # Distortion
    U = runif(n)
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    # Residual estimation
    e.X.est = e.Y.est = vector()
    for (i in 1:n) {
      e.X.est[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
      e.Y.est[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
    }
    
    # Correlation estimate
    cov.e = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
    var.e.X = mean(e.X.est^2) - mean(e.X.est)^2
    var.e.Y = mean(e.Y.est^2) - mean(e.Y.est)^2
    rho.e = cov.e / sqrt(var.e.X * var.e.Y)
    
    # Jackknife pseudo-values
    rho.j = vector()
    for (j in 1:n) {
      jack.X = jack.Y = vector()
      xy.jack = xy.obs[-j, ]
      U.jack = U[-j]
      for (k in 1:(n - 1)) {
        jack.X[k] = xy.jack[k, 1] - log(get.NWK(xy.jack[, 1], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 1])))
        jack.Y[k] = xy.jack[k, 2] - log(get.NWK(xy.jack[, 2], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 2])))
      }
      cov.jack = mean(jack.X * jack.Y) - mean(jack.X) * mean(jack.Y)
      var.X.jack = mean(jack.X^2) - mean(jack.X)^2
      var.Y.jack = mean(jack.Y^2) - mean(jack.Y)^2
      rho.j[j] = n * rho.e - (n - 1) * (cov.jack / sqrt(var.X.jack * var.Y.jack))
    }
    
    # TJEL transformation and test
    trans.rho.j = tanh(rho.j)
    wni = trans.rho.j - mean(trans.rho.j)
    test.tjel = el.test(wni, 0)
    llr = test.tjel$'-2LLR'
    
    # Confidence interval
    ci = findci(trans.rho.j)
    lower.tjel[jj] = tanh(ci$Low)
    upper.tjel[jj] = tanh(ci$Up)
    length.tjel[jj] = upper.tjel[jj] - lower.tjel[jj]
    coverage.tjel[jj] = llr < 1.96^2
  }
  
  results[ii, 1] = mean(lower.tjel)
  results[ii, 2] = mean(upper.tjel)
  results[ii, 3] = mean(length.tjel)
  results[ii, 4] = mean(coverage.tjel)
}

results

