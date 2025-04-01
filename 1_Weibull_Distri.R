# Example 1 AJEL - Weibull Distribution for U
library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 75
rho = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.6)
iter = 100

######## Define Matrix of Results ###################

coverage.jel = vector()
length.jel = vector()
upper.jel = vector()
lower.jel = vector()
results = matrix(NA, nrow = length(rho), ncol = 4)

########## FUNCTIONS #########################

get.NWK = function(x, u, small.u) {
  bw = density(u, kernel = c("epanechnikov"))$bw
  KX = vector()
  K = vector()
  for (j in 1:length(u)) {
    d = (u[j] - small.u) / bw
    KX[j] = exp(x[j]) * (bw^(-1) * max(0, 0.75 * (1 - d^2)))
    K[j] = (bw^(-1) * max(0, 0.75 * (1 - d^2)))
  }
  result = sum(KX) / sum(K)
  return(result)
}

findci <- function(x.vector, AJEL = FALSE) {
  cifunc <- function(ci.val, x.vector, AJEL = AJEL) {
    if (AJEL == TRUE) {
      x.vector = x.vector - ci.val
      x.vector[length(x.vector) + 1] = -0.5 * log(length(x.vector)) * mean(x.vector)
      el.test(x.vector, 0)
    } else {
      el.test(x.vector - ci.val, 0)
    }
  }
  
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector, AJEL = AJEL)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector, AJEL = AJEL)$Up
  return(list(Up = U, Low = L))
}

############ ITERATIONS #######################

for (ii in 1:length(rho)) {
  
  for (jj in 1:iter) {
    
    ####### DATA GENERATION #####################
    
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(-1, 1), Sigma = sigma)
    colnames(data) <- c("X", "Y")
    
    U = rweibull(n, shape = 1.2, scale = 1)  # Weibull distribution for U
    
    ######## DISTORTING FUNCTION #################
    
    phi_X = U^2 - 1.354  # E[U^2] ≈ 1.354, so E[phi_X] = 0
    phi_Y = cos(U) - 0.294  # E[cos(U)] ≈ 0.294, so E[phi_Y] = 0
    
    ###### OBSERVED DATA WITH MULTIPLICATIVE DISTORTION ###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    ####### CALCULATING ESTIMATORS ##############
    
    # RESIDUAL BASED ESTIMATOR ###
    
    e.Y.est = vector()
    for (i in 1:n) {
      e.Y.est[i] = xy.obs[, 2][i] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
    }
    
    e.X.est = vector()
    for (i in 1:n) {
      e.X.est[i] = xy.obs[, 1][i] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
    }
    
    cov.e.est = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
    sig.e.X = mean(e.X.est * e.X.est) - mean(e.X.est) * mean(e.X.est)
    sig.e.Y = mean(e.Y.est * e.Y.est) - mean(e.Y.est) * mean(e.Y.est)
    
    rho.e.est = (cov.e.est) / (sqrt(sig.e.X * sig.e.Y))
    
    ############ JACKKNIFING #######################
    
    rho.j = vector()
    
    for (j in 1:n) {
      
      xy.jack = xy.obs[-j, ]
      U.jack = U[-j]
      
      Eix.j = vector()
      Eiy.j = vector()
      
      for (k in 1:(n-1)) {
        
        Eix.j[k] = xy.jack[, 1][k] - log(get.NWK(xy.jack[, 1], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 1])))
        Eiy.j[k] = xy.jack[, 2][k] - log(get.NWK(xy.jack[, 2], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 2])))
      }
      
      cov.e.est.jack = mean(Eix.j * Eiy.j) - mean(Eix.j) * mean(Eiy.j)
      sig.e.X.jack = mean(Eix.j * Eix.j) - mean(Eix.j) * mean(Eix.j)
      sig.e.Y.jack = mean(Eiy.j * Eiy.j) - mean(Eiy.j) * mean(Eiy.j)
      
      rho.e.est.jack = (cov.e.est.jack) / (sqrt(sig.e.X.jack * sig.e.Y.jack))
      
      rho.j[j] = n * rho.e.est - (n-1) * rho.e.est.jack
    }
    
    ######### JACKKNIFE EMPIRICAL LIKELIHOOD #################
    
    wni = rho.j - rho[ii]
    wni.ci = rho.j - mean(rho.j)
    
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




# Example 1 EL - Weibull Distribution for U
library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 75
rho = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.6)
iter = 100

######## Define Matrix of Results ###################

coverage.el = vector()
length.el = vector()
upper.el = vector()
lower.el = vector()
results = matrix(NA, nrow = length(rho), ncol = 4)

########## FUNCTIONS #########################

get.NWK = function(x, u, small.u) {
  bw = density(u, kernel = c("epanechnikov"))$bw
  KX = vector()
  K = vector()
  for (j in 1:length(u)) {
    d = (u[j] - small.u) / bw
    KX[j] = exp(x[j]) * (bw^(-1) * max(0, 0.75 * (1 - d^2)))
    K[j] = (bw^(-1) * max(0, 0.75 * (1 - d^2)))
  }
  result = sum(KX) / sum(K)
  return(result)
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    el.test(x.vector - ci.val, 0)
  }
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up = U, Low = L))
}

############ ITERATIONS #######################

for (ii in 1:length(rho)) {
  
  for (jj in 1:iter) {
    
    ####### DATA GENERATION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(-1, 1), Sigma = sigma)
    colnames(data) <- c("X", "Y")
    
    U = rweibull(n, shape = 1.2, scale = 1)
    
    ######## DISTORTION FUNCTION #################
    phi_X = U^2 - 1.354    # E[U^2] ≈ 1.354
    phi_Y = cos(U) - 0.294 # E[cos(U)] ≈ 0.294
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    ####### RESIDUALS ##################
    e.X = e.Y = numeric(n)
    for (i in 1:n) {
      e.X[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
      e.Y[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
    }
    
    ######### EMPIRICAL LIKELIHOOD #################
    estimating.eq = (e.X - mean(e.X)) * (e.Y - mean(e.Y))
    
    il.eltest = el.test(estimating.eq, 0)
    il.el = il.eltest$'-2LLR'
    
    ci.el = findci(estimating.eq)
    upper.el[jj] = ci.el$Up
    lower.el[jj] = ci.el$Low
    length.el[jj] = upper.el[jj] - lower.el[jj]
    coverage.el[jj] = il.el < 1.96^2
  }
  
  results[ii, 1] = mean(lower.el)
  results[ii, 2] = mean(upper.el)
  results[ii, 3] = mean(length.el)
  results[ii, 4] = sum(coverage.el) / iter
}

results





# Example 1 JEL - Weibull Distribution for U
library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 75
rho = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.6)
iter = 100

######## Define Matrix of Results ###################

coverage.jel = vector()
length.jel = vector()
upper.jel = vector()
lower.jel = vector()
results = matrix(NA, nrow = length(rho), ncol = 4)

########## FUNCTIONS #########################

get.NWK = function(x, u, small.u) {
  bw = density(u, kernel = c("epanechnikov"))$bw
  KX = vector()
  K = vector()
  for (j in 1:length(u)) {
    d = (u[j] - small.u) / bw
    KX[j] = exp(x[j]) * (bw^(-1) * max(0, 0.75 * (1 - d^2)))
    K[j] = (bw^(-1) * max(0, 0.75 * (1 - d^2)))
  }
  result = sum(KX) / sum(K)
  return(result)
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    el.test(x.vector - ci.val, 0)
  }
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up = U, Low = L))
}

############ ITERATIONS #######################

for (ii in 1:length(rho)) {
  
  for (jj in 1:iter) {
    
    ####### DATA GENERATION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(-1, 1), Sigma = sigma)
    colnames(data) <- c("X", "Y")
    
    U = rweibull(n, shape = 1.2, scale = 1)
    
    ######## DISTORTION FUNCTION #################
    phi_X = U^2 - 1.354    # Centered distortion
    phi_Y = cos(U) - 0.294 # Centered distortion
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    ####### RESIDUALS ##################
    e.X = e.Y = numeric(n)
    for (i in 1:n) {
      e.X[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
      e.Y[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
    }
    
    cov.e = mean(e.X * e.Y) - mean(e.X) * mean(e.Y)
    var.e.X = mean(e.X^2) - mean(e.X)^2
    var.e.Y = mean(e.Y^2) - mean(e.Y)^2
    rho.hat = cov.e / sqrt(var.e.X * var.e.Y)
    
    ########## JACKKNIFE PSEUDO-VALUES ########################
    V = numeric(n)
    for (j in 1:n) {
      e.X.j = e.X[-j]
      e.Y.j = e.Y[-j]
      cov.j = mean(e.X.j * e.Y.j) - mean(e.X.j) * mean(e.Y.j)
      var.X.j = mean(e.X.j^2) - mean(e.X.j)^2
      var.Y.j = mean(e.Y.j^2) - mean(e.Y.j)^2
      rho.j = cov.j / sqrt(var.X.j * var.Y.j)
      V[j] = n * rho.hat - (n - 1) * rho.j
    }
    
    ########## JEL ESTIMATION ##################
    W = V - rho[ii]
    il.jel = el.test(W, 0)$'-2LLR'
    
    ci.jel = findci(V)
    upper.jel[jj] = ci.jel$Up
    lower.jel[jj] = ci.jel$Low
    length.jel[jj] = upper.jel[jj] - lower.jel[jj]
    coverage.jel[jj] = il.jel < 1.96^2
  }
  
  results[ii, 1] = mean(lower.jel)
  results[ii, 2] = mean(upper.jel)
  results[ii, 3] = mean(length.jel)
  results[ii, 4] = sum(coverage.jel) / iter
}

results




# Example 1 JEL - Weibull Distribution for U
library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 75
rho = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.6)
iter = 100

######## Define Matrix of Results ###################

coverage.jel = vector()
length.jel = vector()
upper.jel = vector()
lower.jel = vector()
results = matrix(NA, nrow = length(rho), ncol = 4)

########## FUNCTIONS #########################

get.NWK = function(x, u, small.u) {
  bw = density(u, kernel = c("epanechnikov"))$bw
  KX = vector()
  K = vector()
  for (j in 1:length(u)) {
    d = (u[j] - small.u) / bw
    KX[j] = exp(x[j]) * (bw^(-1) * max(0, 0.75 * (1 - d^2)))
    K[j] = (bw^(-1) * max(0, 0.75 * (1 - d^2)))
  }
  result = sum(KX) / sum(K)
  return(result)
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    el.test(x.vector - ci.val, 0)
  }
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up = U, Low = L))
}

############ ITERATIONS #######################

for (ii in 1:length(rho)) {
  
  for (jj in 1:iter) {
    
    ####### DATA GENERATION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(-1, 1), Sigma = sigma)
    colnames(data) <- c("X", "Y")
    
    U = rweibull(n, shape = 1.2, scale = 1)
    
    ######## DISTORTION FUNCTION #################
    phi_X = U^2 - 1.354    # Centered distortion
    phi_Y = cos(U) - 0.294 # Centered distortion
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    ####### RESIDUALS ##################
    e.X = e.Y = numeric(n)
    for (i in 1:n) {
      e.X[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
      e.Y[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
    }
    
    cov.e = mean(e.X * e.Y) - mean(e.X) * mean(e.Y)
    var.e.X = mean(e.X^2) - mean(e.X)^2
    var.e.Y = mean(e.Y^2) - mean(e.Y)^2
    rho.hat = cov.e / sqrt(var.e.X * var.e.Y)
    
    ########## JACKKNIFE PSEUDO-VALUES ########################
    V = numeric(n)
    for (j in 1:n) {
      e.X.j = e.X[-j]
      e.Y.j = e.Y[-j]
      cov.j = mean(e.X.j * e.Y.j) - mean(e.X.j) * mean(e.Y.j)
      var.X.j = mean(e.X.j^2) - mean(e.X.j)^2
      var.Y.j = mean(e.Y.j^2) - mean(e.Y.j)^2
      rho.j = cov.j / sqrt(var.X.j * var.Y.j)
      V[j] = n * rho.hat - (n - 1) * rho.j
    }
    
    ########## JEL ESTIMATION ##################
    W = V - rho[ii]
    il.jel = el.test(W, 0)$'-2LLR'
    
    ci.jel = findci(V)
    upper.jel[jj] = ci.jel$Up
    lower.jel[jj] = ci.jel$Low
    length.jel[jj] = upper.jel[jj] - lower.jel[jj]
    coverage.jel[jj] = il.jel < 1.96^2
  }
  
  results[ii, 1] = mean(lower.jel)
  results[ii, 2] = mean(upper.jel)
  results[ii, 3] = mean(length.jel)
  results[ii, 4] = sum(coverage.jel) / iter
}

results





# Example 1 MJEL - Weibull Distribution for U
library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 75
rho = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.6)
iter = 100

######## Define Matrix of Results ###################

coverage.mjel = vector()
length.mjel = vector()
upper.mjel = vector()
lower.mjel = vector()
results = matrix(NA, nrow = length(rho), ncol = 4)

########## FUNCTIONS #########################

get.NWK = function(x, u, small.u) {
  bw = density(u, kernel = c("epanechnikov"))$bw
  KX = vector()
  K = vector()
  for (j in 1:length(u)) {
    d = (u[j] - small.u) / bw
    KX[j] = exp(x[j]) * (bw^(-1) * max(0, 0.75 * (1 - d^2)))
    K[j] = (bw^(-1) * max(0, 0.75 * (1 - d^2)))
  }
  result = sum(KX) / sum(K)
  return(result)
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    el.test(x.vector - ci.val, 0)
  }
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up = U, Low = L))
}

############ ITERATIONS #######################

for (ii in 1:length(rho)) {
  
  for (jj in 1:iter) {
    
    ####### DATA GENERATION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(-1, 1), Sigma = sigma)
    colnames(data) <- c("X", "Y")
    
    U = rweibull(n, shape = 1.2, scale = 1)
    
    ######## DISTORTION FUNCTION #################
    phi_X = U^2 - 1.354    # E[U^2] ≈ 1.354
    phi_Y = cos(U) - 0.294 # E[cos(U)] ≈ 0.294
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    ####### RESIDUALS ##################
    e.X = e.Y = numeric(n)
    for (i in 1:n) {
      e.X[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
      e.Y[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
    }
    
    cov.e = mean(e.X * e.Y) - mean(e.X) * mean(e.Y)
    var.e.X = mean(e.X^2) - mean(e.X)^2
    var.e.Y = mean(e.Y^2) - mean(e.Y)^2
    rho.hat = cov.e / sqrt(var.e.X * var.e.Y)
    
    ########## JACKKNIFE PSEUDO-VALUES ########################
    V = numeric(n)
    for (j in 1:n) {
      e.X.j = e.X[-j]
      e.Y.j = e.Y[-j]
      cov.j = mean(e.X.j * e.Y.j) - mean(e.X.j) * mean(e.Y.j)
      var.X.j = mean(e.X.j^2) - mean(e.X.j)^2
      var.Y.j = mean(e.Y.j^2) - mean(e.Y.j)^2
      rho.j = cov.j / sqrt(var.X.j * var.Y.j)
      V[j] = n * rho.hat - (n - 1) * rho.j
    }
    
    ########## MJEL ESTIMATION ##################
    W.centered = V - mean(V)
    il.mjel = el.test(W.centered, 0)$'-2LLR'
    
    ci.mjel = findci(V)
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



# Example 1 TJEL - Weibull Distribution for U
library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 75
rho = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.6)
iter = 100

######## Define Matrix of Results ###################

coverage.tjel = vector()
length.tjel = vector()
upper.tjel = vector()
lower.tjel = vector()
results = matrix(NA, nrow = length(rho), ncol = 4)

########## FUNCTIONS #########################

get.NWK = function(x, u, small.u) {
  bw = density(u, kernel = c("epanechnikov"))$bw
  KX = vector()
  K = vector()
  for (j in 1:length(u)) {
    d = (u[j] - small.u) / bw
    KX[j] = exp(x[j]) * (bw^(-1) * max(0, 0.75 * (1 - d^2)))
    K[j] = (bw^(-1) * max(0, 0.75 * (1 - d^2)))
  }
  result = sum(KX) / sum(K)
  return(result)
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    el.test(x.vector - ci.val, 0)
  }
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up = U, Low = L))
}

############ ITERATIONS #######################

for (ii in 1:length(rho)) {
  
  for (jj in 1:iter) {
    
    ####### DATA GENERATION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(-1, 1), Sigma = sigma)
    colnames(data) <- c("X", "Y")
    
    U = rweibull(n, shape = 1.2, scale = 1)
    
    ######## DISTORTION FUNCTION #################
    phi_X = U^2 - 1.354
    phi_Y = cos(U) - 0.294
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    ####### RESIDUALS ##################
    e.X = e.Y = numeric(n)
    for (i in 1:n) {
      e.X[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
      e.Y[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
    }
    
    cov.e = mean(e.X * e.Y) - mean(e.X) * mean(e.Y)
    var.e.X = mean(e.X^2) - mean(e.X)^2
    var.e.Y = mean(e.Y^2) - mean(e.Y)^2
    rho.hat = cov.e / sqrt(var.e.X * var.e.Y)
    
    ########## JACKKNIFE PSEUDO-VALUES ########################
    V = numeric(n)
    for (j in 1:n) {
      e.X.j = e.X[-j]
      e.Y.j = e.Y[-j]
      cov.j = mean(e.X.j * e.Y.j) - mean(e.X.j) * mean(e.Y.j)
      var.X.j = mean(e.X.j^2) - mean(e.X.j)^2
      var.Y.j = mean(e.Y.j^2) - mean(e.Y.j)^2
      rho.j = cov.j / sqrt(var.X.j * var.Y.j)
      V[j] = n * rho.hat - (n - 1) * rho.j
    }
    
    ########## TRANSFORMATION: FISHER'S Z #####################
    V[V <= -0.9999] <- -0.9999
    V[V >=  0.9999] <-  0.9999
    V.transformed <- atanh(V)
    
    ########## TJEL ESTIMATION ##################
    W.centered = V.transformed - mean(V.transformed)
    il.tjel = el.test(W.centered, 0)$'-2LLR'
    
    ci.tjel = findci(V.transformed)
    upper.tjel[jj] = tanh(ci.tjel$Up)
    lower.tjel[jj] = tanh(ci.tjel$Low)
    length.tjel[jj] = upper.tjel[jj] - lower.tjel[jj]
    coverage.tjel[jj] = il.tjel < 1.96^2
  }
  
  results[ii, 1] = mean(lower.tjel)
  results[ii, 2] = mean(upper.tjel)
  results[ii, 3] = mean(length.tjel)
  results[ii, 4] = sum(coverage.tjel) / iter
}

results

