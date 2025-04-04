install.packages("emplik")
install.packages("kedd")

library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 50
rho = c(-0.9, -0.5, 0, 0.5, 0.9)
iter = 50

######## Define Matrix of Results ###########

coverage.jel = vector()
length.jel = vector()
upper.jel = vector()
lower.jel = vector()
results = matrix(NA, nrow = 5, ncol = 4)

########## FUNCTIONS #########################

get.NWK = function(x, u, small.u){
  bw = density(u, kernel = c("epanechnikov"))$bw * 1.2
  KX = vector()
  K = vector()
  for(j in 1:length(u)){
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

############ ITERATIONS ######################

for (ii in 1:5){
  for (jj in 1:iter){
    
    ####### DATA GENERATION ####################
    
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma) 
    colnames(data) <- c("X", "Y")
    
    U = runif(n, -1, 1)
    
    ######## MULTIPLICATIVE DISTORTION FUNCTION ###############
    
    phi_X = 1 + 0.25 * sin(2 * pi * U)  # Range ~ [0.75, 1.25]
    phi_Y = 1 + 0.15 * (-U + log(2 / (exp(1) - exp(-1))))  # Range ~ [~0.85, ~1.15]
    
    ###### OBSERVED DATA ########################
    
    xy.obs <- data * cbind(phi_X, phi_Y)
    
    ####### CALCULATING RESIDUAL-BASED ESTIMATORS ###########
    
    e.Y.est = vector()
    for (i in 1:n){
      e.Y.est[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
    }
    
    e.X.est = vector()
    for (i in 1:n){
      e.X.est[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
    }
    
    cov.e.est = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
    sig.e.X = mean(e.X.est^2) - mean(e.X.est)^2
    sig.e.Y = mean(e.Y.est^2) - mean(e.Y.est)^2
    
    rho.e.est = cov.e.est / sqrt(sig.e.X * sig.e.Y)
    
    ############ JACKKNIFING #####################
    
    rho.j = vector()
    
    for (j in 1:n){
      xy.jack = xy.obs[-j, ]
      U.jack = U[-j]
      
      Eix.j = vector()
      Eiy.j = vector()
      
      for (k in 1:(n - 1)){
        Eix.j[k] = xy.jack[k, 1] - log(get.NWK(xy.jack[, 1], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 1])))
        Eiy.j[k] = xy.jack[k, 2] - log(get.NWK(xy.jack[, 2], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 2])))
      }
      
      cov.e.est.jack = mean(Eix.j * Eiy.j) - mean(Eix.j) * mean(Eiy.j)
      sig.e.X.jack = mean(Eix.j^2) - mean(Eix.j)^2
      sig.e.Y.jack = mean(Eiy.j^2) - mean(Eiy.j)^2
      
      rho.e.est.jack = cov.e.est.jack / sqrt(sig.e.X.jack * sig.e.Y.jack)
      
      rho.j[j] = n * rho.e.est - (n - 1) * rho.e.est.jack
    }
    
    ######### JACKKNIFE EMPIRICAL LIKELIHOOD #########
    
    wni = rho.j - rho[ii]
    
    il.jel = el.test(wni, 0)$'-2LLR'
    ci.jel = findci(rho.j)
    
    upper.jel[jj] = ci.jel$Up
    lower.jel[jj] = ci.jel$Low
    length.jel[jj] = upper.jel[jj] - lower.jel[jj]
    coverage.jel[jj] = il.jel < qchisq(0.975, 1)
  }
  
  results[ii, 1] = mean(lower.jel)
  results[ii, 2] = mean(upper.jel)
  results[ii, 3] = mean(length.jel)
  results[ii, 4] = sum(coverage.jel) / iter
}

results
