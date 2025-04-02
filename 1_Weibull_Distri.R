# #example 1 AJEL with Weibull-distributed data

library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 75
rho = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.6)
iter = 100

########Define Matrix of Results ###################

coverage.jel = vector()
length.jel = vector()
upper.jel = vector()
lower.jel = vector()
results = matrix(NA, nrow = length(rho), ncol = 4)

########## FUNCTIONS #########################

get.NWK = function(x, u, small.u) {
  bw = density(U, kernel = c("epanechnikov"))$bw
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
      x.vector[length(x.vector)+1] = -0.5 * log(length(x.vector)) * mean(x.vector)
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
    
    ####### DATA GENERATION using WEIBULL DISTRIBUTION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    norm_data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
    
    U1 <- pnorm(norm_data[, 1])
    U2 <- pnorm(norm_data[, 2])
    
    # Transform to Weibull(2,1)
    X <- qweibull(U1, shape = 2, scale = 1)
    Y <- qweibull(U2, shape = 2, scale = 1)
    
    data <- cbind(X, Y)
    colnames(data) <- c("X", "Y")
    
    U = runif(n, 0, 1)
    
    ######## DISTORTING FUNCTION #################
    
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    
    ###### OBSERVED DATA WITH MULTIPLICATIVE DISTORTION ###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    ####### CALCULATING ESTIMATORS ####################
    
    e.Y.est = vector()
    e.X.est = vector()
    for (i in 1:n) {
      e.Y.est[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
      e.X.est[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
    }
    
    cov.e.est = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
    sig.e.X = mean(e.X.est^2) - mean(e.X.est)^2
    sig.e.Y = mean(e.Y.est^2) - mean(e.Y.est)^2
    
    rho.e.est = cov.e.est / sqrt(sig.e.X * sig.e.Y)
    
    ############ JACKKNIFING #######################
    
    rho.j = vector()
    
    for (j in 1:n) {
      
      xy.jack = xy.obs[-j, ]
      U.jack = U[-j]
      
      Eix.j = vector()
      Eiy.j = vector()
      
      for (k in 1:(n - 1)) {
        Eix.j[k] = xy.jack[k, 1] - log(get.NWK(xy.jack[, 1], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 1])))
        Eiy.j[k] = xy.jack[k, 2] - log(get.NWK(xy.jack[, 2], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 2])))
      }
      
      cov.e.est.jack = mean(Eix.j * Eiy.j) - mean(Eix.j) * mean(Eiy.j)
      sig.e.X.jack = mean(Eix.j^2) - mean(Eix.j)^2
      sig.e.Y.jack = mean(Eiy.j^2) - mean(Eiy.j)^2
      
      rho.e.est.jack = cov.e.est.jack / sqrt(sig.e.X.jack * sig.e.Y.jack)
      
      rho.j[j] = n * rho.e.est - (n - 1) * rho.e.est.jack
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


# EL method with Weibull-distributed data

library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 75
rho = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.6)
iter = 100

coverage.el = vector()
length.el = vector()
upper.el = vector()
lower.el = vector()
results = matrix(NA, nrow = length(rho), ncol = 4)

########## FUNCTIONS #########################

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
    
    ####### DATA GENERATION using WEIBULL DISTRIBUTION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    norm_data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
    
    U1 <- pnorm(norm_data[, 1])
    U2 <- pnorm(norm_data[, 2])
    
    # Transform to Weibull(2,1)
    X <- qweibull(U1, shape = 2, scale = 1)
    Y <- qweibull(U2, shape = 2, scale = 1)
    
    data <- cbind(X, Y)
    colnames(data) <- c("X", "Y")
    
    U = runif(n, 0, 1)
    
    ######## DISTORTING FUNCTION #################
    
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    
    ###### OBSERVED DATA WITH MULTIPLICATIVE DISTORTION ###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    ####### RESIDUAL-BASED ESTIMATOR ####################
    
    e.X.est = e.Y.est = vector()
    for (i in 1:n) {
      e.X.est[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
      e.Y.est[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
    }
    
    # Compute pseudo-values for EL test
    cov.e = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
    sig.e.X = mean(e.X.est^2) - mean(e.X.est)^2
    sig.e.Y = mean(e.Y.est^2) - mean(e.Y.est)^2
    rho.e = cov.e / sqrt(sig.e.X * sig.e.Y)
    
    wni = (e.X.est - mean(e.X.est)) * (e.Y.est - mean(e.Y.est)) / sqrt(sig.e.X * sig.e.Y)
    wni = wni - rho[ii]
    
    il.eltest = el.test(wni, 0)
    llr = il.eltest$'-2LLR'
    
    ci.el = findci(wni)
    upper.el[jj] = ci.el$Up
    lower.el[jj] = ci.el$Low
    length.el[jj] = upper.el[jj] - lower.el[jj]
    coverage.el[jj] = llr < 1.96^2
  }
  
  results[ii, 1] = mean(lower.el)
  results[ii, 2] = mean(upper.el)
  results[ii, 3] = mean(length.el)
  results[ii, 4] = sum(coverage.el) / iter
}

results




# JEL method with Weibull-distributed data

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
  KX = K = vector()
  for (j in 1:length(u)) {
    d = (u[j] - small.u) / bw
    KX[j] = exp(x[j]) * (bw^(-1) * max(0, 0.75 * (1 - d^2)))
    K[j] = (bw^(-1) * max(0, 0.75 * (1 - d^2)))
  }
  return(sum(KX) / sum(K))
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
    
    ####### DATA GENERATION using WEIBULL DISTRIBUTION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    norm_data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
    
    U1 <- pnorm(norm_data[, 1])
    U2 <- pnorm(norm_data[, 2])
    
    X <- qweibull(U1, shape = 2, scale = 1)
    Y <- qweibull(U2, shape = 2, scale = 1)
    
    data <- cbind(X, Y)
    colnames(data) <- c("X", "Y")
    
    U = runif(n, 0, 1)
    
    ######## DISTORTING FUNCTION #################
    
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    
    ###### OBSERVED DATA WITH MULTIPLICATIVE DISTORTION ###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    ####### CALCULATING RESIDUAL-BASED ESTIMATORS ####################
    
    e.X.est = e.Y.est = vector()
    for (i in 1:n) {
      e.X.est[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
      e.Y.est[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
    }
    
    cov.e = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
    sig.e.X = mean(e.X.est^2) - mean(e.X.est)^2
    sig.e.Y = mean(e.Y.est^2) - mean(e.Y.est)^2
    rho.e = cov.e / sqrt(sig.e.X * sig.e.Y)
    
    ########## JACKKNIFING ###########################
    
    rho.j = vector()
    for (j in 1:n) {
      xy.jack = xy.obs[-j, ]
      U.jack = U[-j]
      
      jack.X = jack.Y = vector()
      for (k in 1:(n - 1)) {
        jack.X[k] = xy.jack[k, 1] - log(get.NWK(xy.jack[, 1], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 1])))
        jack.Y[k] = xy.jack[k, 2] - log(get.NWK(xy.jack[, 2], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 2])))
      }
      
      cov.j = mean(jack.X * jack.Y) - mean(jack.X) * mean(jack.Y)
      var.X.j = mean(jack.X^2) - mean(jack.X)^2
      var.Y.j = mean(jack.Y^2) - mean(jack.Y)^2
      
      rho.j[j] = n * rho.e - (n - 1) * (cov.j / sqrt(var.X.j * var.Y.j))
    }
    
    ########## JEL EMPIRICAL LIKELIHOOD ##################
    
    wni = rho.j - rho[ii]
    il.jel = el.test(wni, 0)$'-2LLR'
    ci.jel = findci(rho.j)
    
    upper.jel[jj] = ci.jel$Up
    lower.jel[jj] = ci.jel$Low
    length.jel[jj] = upper.jel[jj] - lower.jel[jj]
    coverage.jel[jj] = il.jel < 1.96^2
  }
  
  results[ii, 1] = mean(lower.jel)
  results[ii, 2] = mean(upper.jel)
  results[ii, 3] = mean(length.jel)
  results[ii, 4] = mean(coverage.jel)
}

results





# MJEL method with Weibull-distributed data

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
  KX = K = vector()
  for (j in 1:length(u)) {
    d = (u[j] - small.u) / bw
    KX[j] = exp(x[j]) * (bw^(-1) * max(0, 0.75 * (1 - d^2)))
    K[j] = (bw^(-1) * max(0, 0.75 * (1 - d^2)))
  }
  return(sum(KX) / sum(K))
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
    
    ####### DATA GENERATION using WEIBULL DISTRIBUTION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    norm_data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
    
    U1 <- pnorm(norm_data[, 1])
    U2 <- pnorm(norm_data[, 2])
    
    X <- qweibull(U1, shape = 2, scale = 1)
    Y <- qweibull(U2, shape = 2, scale = 1)
    
    data <- cbind(X, Y)
    colnames(data) <- c("X", "Y")
    
    U = runif(n, 0, 1)
    
    ######## DISTORTING FUNCTION #################
    
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    
    ###### OBSERVED DATA WITH MULTIPLICATIVE DISTORTION ###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    ####### CALCULATING RESIDUAL-BASED ESTIMATORS ####################
    
    e.X.est = e.Y.est = vector()
    for (i in 1:n) {
      e.X.est[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
      e.Y.est[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
    }
    
    cov.e = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
    sig.e.X = mean(e.X.est^2) - mean(e.X.est)^2
    sig.e.Y = mean(e.Y.est^2) - mean(e.Y.est)^2
    rho.e = cov.e / sqrt(sig.e.X * sig.e.Y)
    
    ########## JACKKNIFING ###########################
    
    rho.j = vector()
    for (j in 1:n) {
      xy.jack = xy.obs[-j, ]
      U.jack = U[-j]
      
      jack.X = jack.Y = vector()
      for (k in 1:(n - 1)) {
        jack.X[k] = xy.jack[k, 1] - log(get.NWK(xy.jack[, 1], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 1])))
        jack.Y[k] = xy.jack[k, 2] - log(get.NWK(xy.jack[, 2], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 2])))
      }
      
      cov.j = mean(jack.X * jack.Y) - mean(jack.X) * mean(jack.Y)
      var.X.j = mean(jack.X^2) - mean(jack.X)^2
      var.Y.j = mean(jack.Y^2) - mean(jack.Y)^2
      
      rho.j[j] = n * rho.e - (n - 1) * (cov.j / sqrt(var.X.j * var.Y.j))
    }
    
    ########## MJEL EMPIRICAL LIKELIHOOD ##################
    
    rho.bar = mean(rho.j)
    wni = rho.j - rho.bar
    
    il.mjel = el.test(wni, 0)$'-2LLR'
    ci.mjel = findci(rho.j)
    
    upper.mjel[jj] = ci.mjel$Up
    lower.mjel[jj] = ci.mjel$Low
    length.mjel[jj] = upper.mjel[jj] - lower.mjel[jj]
    coverage.mjel[jj] = il.mjel < 1.96^2
  }
  
  results[ii, 1] = mean(lower.mjel)
  results[ii, 2] = mean(upper.mjel)
  results[ii, 3] = mean(length.mjel)
  results[ii, 4] = mean(coverage.mjel)
}

results



# TJEL method with Weibull-distributed data

library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 75
rho = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.6)
iter = 100

coverage.tjel = vector()
length.tjel = vector()
upper.tjel = vector()
lower.tjel = vector()
results = matrix(NA, nrow = length(rho), ncol = 4)

########## FUNCTIONS #########################

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

findci_transformed <- function(x.vector) {
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
    
    ####### DATA GENERATION using WEIBULL DISTRIBUTION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    norm_data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
    
    U1 <- pnorm(norm_data[, 1])
    U2 <- pnorm(norm_data[, 2])
    
    X <- qweibull(U1, shape = 2, scale = 1)
    Y <- qweibull(U2, shape = 2, scale = 1)
    
    data <- cbind(X, Y)
    colnames(data) <- c("X", "Y")
    
    U = runif(n, 0, 1)
    
    ######## DISTORTING FUNCTION #################
    
    phi_X = 1.25 - 3 * (U - 0.5)^2
    phi_Y = 1 + 0.5 * cos(2 * pi * U)
    
    ###### OBSERVED DATA WITH MULTIPLICATIVE DISTORTION ###########
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    ####### CALCULATING RESIDUAL-BASED ESTIMATORS ####################
    
    e.X.est = e.Y.est = vector()
    for (i in 1:n) {
      e.X.est[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
      e.Y.est[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
    }
    
    cov.e = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
    sig.e.X = mean(e.X.est^2) - mean(e.X.est)^2
    sig.e.Y = mean(e.Y.est^2) - mean(e.Y.est)^2
    rho.e = cov.e / sqrt(sig.e.X * sig.e.Y)
    
    ########## JACKKNIFING ###########################
    
    rho.j = vector()
    for (j in 1:n) {
      xy.jack = xy.obs[-j, ]
      U.jack = U[-j]
      
      jack.X = jack.Y = vector()
      for (k in 1:(n - 1)) {
        jack.X[k] = xy.jack[k, 1] - log(get.NWK(xy.jack[, 1], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 1])))
        jack.Y[k] = xy.jack[k, 2] - log(get.NWK(xy.jack[, 2], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 2])))
      }
      
      cov.j = mean(jack.X * jack.Y) - mean(jack.X) * mean(jack.Y)
      var.X.j = mean(jack.X^2) - mean(jack.X)^2
      var.Y.j = mean(jack.Y^2) - mean(jack.Y)^2
      
      rho.j[j] = n * rho.e - (n - 1) * (cov.j / sqrt(var.X.j * var.Y.j))
    }
    
    ########## TJEL EMPIRICAL LIKELIHOOD ##################
    
    # Transform pseudo-values using tanh
    trans_rho_j = tanh(rho.j)
    wni = trans_rho_j - mean(trans_rho_j)
    
    il.tjel = el.test(wni, 0)$'-2LLR'
    ci.tjel = findci_transformed(trans_rho_j)
    
    # Back-transform CI endpoints
    lower.tjel[jj] = tanh(ci.tjel$Low)
    upper.tjel[jj] = tanh(ci.tjel$Up)
    length.tjel[jj] = upper.tjel[jj] - lower.tjel[jj]
    coverage.tjel[jj] = il.tjel < 1.96^2
  }
  
  results[ii, 1] = mean(lower.tjel)
  results[ii, 2] = mean(upper.tjel)
  results[ii, 3] = mean(length.tjel)
  results[ii, 4] = mean(coverage.tjel)
}

results



