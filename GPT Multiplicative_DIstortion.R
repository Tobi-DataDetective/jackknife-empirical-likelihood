library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 50
rho = c(-0.9, -0.5, 0, 0.5, 0.9)
iter = 2000

######## Define Matrix of Results ###################

coverage.jel = vector()
length.jel = vector()
upper.jel = vector()
lower.jel = vector()
results = matrix(NA, nrow = 5, ncol = 4)

########## FUNCTIONS #########################

get.NWK = function(x, u, small.u){
  bw = density(U, kernel = "epanechnikov")$bw
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

for (ii in 1:5){
  
  for (jj in 1:iter){
    
    ####### DATA GENERATION #####################
    
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
    colnames(data) <- c("X", "Y")
    
    U = runif(n, -1, 1)  # Multiplicative distortion source
    
    ######## DISTORTING FUNCTION #################
    
    psi_U = exp(0.5 * sin(pi * U))    # E[psi(U)] ≈ 1
    phi_U = exp(0.5 * cos(pi * U))    # E[phi(U)] ≈ 1
    
    ####### OBSERVED DATA ########################
    
    xy.obs = data * cbind(psi_U, phi_U)
    
    ####### CALCULATING ESTIMATORS ###############
    
    e.X.est = vector()
    e.Y.est = vector()
    
    for (i in 1:n){
      e.X.est[i] = log(xy.obs[,1][i]) - log(get.NWK(log(xy.obs[,1]), U, U[i])) + log(mean(exp(log(xy.obs[,1]))))
      e.Y.est[i] = log(xy.obs[,2][i]) - log(get.NWK(log(xy.obs[,2]), U, U[i])) + log(mean(exp(log(xy.obs[,2]))))
    }
    
    cov.e.est = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
    sig.e.X = var(e.X.est)
    sig.e.Y = var(e.Y.est)
    
    rho.e.est = cov.e.est / sqrt(sig.e.X * sig.e.Y)
    
    ########### JACKKNIFING ######################
    
    rho.j = vector()
    
    for (j in 1:n){
      xy.jack = xy.obs[-j,]
      U.jack = U[-j]
      
      Eix.j = vector()
      Eiy.j = vector()
      
      for (k in 1:(n - 1)){
        Eix.j[k] = log(xy.jack[,1][k]) - log(get.NWK(log(xy.jack[,1]), U.jack, U.jack[k])) + log(mean(exp(log(xy.jack[,1]))))
        Eiy.j[k] = log(xy.jack[,2][k]) - log(get.NWK(log(xy.jack[,2]), U.jack, U.jack[k])) + log(mean(exp(log(xy.jack[,2]))))
      }
      
      cov.e.est.jack = mean(Eix.j * Eiy.j) - mean(Eix.j) * mean(Eiy.j)
      sig.e.X.jack = var(Eix.j)
      sig.e.Y.jack = var(Eiy.j)
      
      rho.e.est.jack = cov.e.est.jack / sqrt(sig.e.X.jack * sig.e.Y.jack)
      
      rho.j[j] = n * rho.e.est - (n - 1) * rho.e.est.jack
    }
    
    ######### JEL CONFIDENCE INTERVAL #################
    
    wni = rho.j - rho[ii]
    
    il.jeltest = el.test(wni, 0)
    il.j = il.jeltest$`-2LLR`
    
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
