# JEL
# install.packages("emplik")
# install.packages("kedd")


library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS########################

n = 50
rho = c(-0.9, -0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7, 0.9)
iter = 50

# Define the distributions for U
distributions = list(
  list(name = "Uniform", gen_U = function(n) runif(n, -1, 1), mean_U = 0, mean_U2 = 1/3, mean_expU = (exp(1) - exp(-1))/2),
  list(name = "Normal", gen_U = function(n) rnorm(n, 0, 1), mean_U = 0, mean_U2 = 1, mean_expU = exp(0.5)),
  list(name = "Beta", gen_U = function(n) rbeta(n, 2, 2) - 0.5, mean_U = 0, mean_U2 = 0.05, mean_expU = (exp(0.5) - exp(-0.5))/2)
)

# List to store results for each distribution and distortion function
all_results = list()

########Define Matrix of Results###################

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

############ITERATIONS OVER DISTRIBUTIONS AND DISTORTION FUNCTIONS#######################

for (distrib_idx in 1:length(distributions)) {
  distrib_name = distributions[[distrib_idx]]$name
  gen_U = distributions[[distrib_idx]]$gen_U
  mean_U = distributions[[distrib_idx]]$mean_U
  mean_U2 = distributions[[distrib_idx]]$mean_U2
  mean_expU = distributions[[distrib_idx]]$mean_expU
  
  # Define distortion functions based on the distribution
  distortion_functions = list(
    list(name = "Advanced_Linear", 
         phi_X = function(U) 0.5 * U + log(2 / (exp(1) - exp(-1))), 
         phi_Y = function(U) -0.5 * U - log(2 / (exp(1) - exp(-1)))),
    list(name = "Advanced_Quadratic", 
         phi_X = function(U) U^2 - mean_U2 + 0.1 * exp(U - mean_expU), 
         phi_Y = function(U) -(U^2 - mean_U2) - 0.1 * exp(U - mean_expU)),
    list(name = "Advanced_Periodic", 
         phi_X = function(U) (1 + 0.5 * U) * sin(pi * U), 
         phi_Y = function(U) -(1 + 0.5 * U) * sin(pi * U))
  )
  
  # List to store results for this distribution
  distrib_results = list()
  
  for (dist_idx in 1:length(distortion_functions)) {
    distortion_name = distortion_functions[[dist_idx]]$name
    phi_X_func = distortion_functions[[dist_idx]]$phi_X
    phi_Y_func = distortion_functions[[dist_idx]]$phi_Y
    
    # Matrix to store results for this distortion function
    results = matrix(NA, nrow = length(rho), ncol = 4)
    colnames(results) = c("Lower", "Upper", "Length", "Coverage")
    rownames(results) = as.character(rho)  # <--- Added this line
    
    
    coverage.jel = vector()
    length.jel = vector()
    upper.jel = vector()
    lower.jel = vector()
    
    for (ii in 1:length(rho)) {
      for (jj in 1:iter) {
        ####### DATA GENERATION#####################
        
        sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
        data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
        colnames(data) <- c("X", "Y")
        
        U = gen_U(n)
        
        ########DISTORTING FUNCTION#################
        
        phi_X = phi_X_func(U)
        phi_Y = phi_Y_func(U)
        
        ######OBSERVED DATA#########################
        
        xy.obs <- data + cbind(phi_X, phi_Y)
        
        #######CALCULATING ESTIMATORS##############
        
        #RESIDUAL BASED ESTIMATOR###
        
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
        
        ############JACKKNIFING#######################
        
        rho.j = vector()
        Eix.j = vector()
        Eiy.j = vector()
        
        for (j in 1:n) {
          xy.jack = xy.obs[-j, ]
          U.jack = U[-j]
          
          Eix.j = vector()
          Eiy.j = vector()
          
          for (k in 1:(n - 1)) {
            Eix.j[k] = xy.jack[, 1][k] - log(get.NWK(xy.jack[, 1], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 1])))
            Eiy.j[k] = xy.jack[, 2][k] - log(get.NWK(xy.jack[, 2], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 2])))
          }
          
          cov.e.est.jack = mean(Eix.j * Eiy.j) - mean(Eix.j) * mean(Eiy.j)
          sig.e.X.jack = mean(Eix.j * Eix.j) - mean(Eix.j) * mean(Eix.j)
          sig.e.Y.jack = mean(Eiy.j * Eiy.j) - mean(Eiy.j) * mean(Eiy.j)
          
          rho.e.est.jack = (cov.e.est.jack) / (sqrt(sig.e.X.jack * sig.e.Y.jack))
          
          rho.j[j] = n * rho.e.est - (n - 1) * rho.e.est.jack
        }
        
        #########JACKKNIFE EMPIRICAL LIKELIHOOD#################
        
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
    
    distrib_results[[distortion_name]] = results
  }
  
  all_results[[distrib_name]] = distrib_results
}

# Print results
all_results





# AJEL
library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 50
rho = c(-0.9, -0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7, 0.9)
iter = 50

distributions = list(
  list(name = "Uniform", gen_U = function(n) runif(n, -1, 1), mean_U = 0, mean_U2 = 1/3, mean_expU = (exp(1) - exp(-1))/2),
  list(name = "Normal", gen_U = function(n) rnorm(n, 0, 1), mean_U = 0, mean_U2 = 1, mean_expU = exp(0.5)),
  list(name = "Beta", gen_U = function(n) rbeta(n, 2, 2) - 0.5, mean_U = 0, mean_U2 = 0.05, mean_expU = (exp(0.5) - exp(-0.5))/2)
)

all_results = list()

get.NWK <- function(x, u, small.u) {
  bw = density(U, kernel = c("epanechnikov"))$bw
  KX = K = numeric(length(u))
  for (j in 1:length(u)) {
    d = (u[j] - small.u) / bw
    weight = bw^(-1) * max(0, 0.75 * (1 - d^2))
    KX[j] = exp(x[j]) * weight
    K[j] = weight
  }
  return(sum(KX) / sum(K))
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    x.adj = x.vector - ci.val
    x.adj[length(x.adj) + 1] = -0.5 * log(length(x.vector)) * mean(x.adj)
    el.test(x.adj, 0)
  }
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up = U, Low = L))
}

for (distrib_idx in seq_along(distributions)) {
  distrib_name = distributions[[distrib_idx]]$name
  gen_U = distributions[[distrib_idx]]$gen_U
  mean_U2 = distributions[[distrib_idx]]$mean_U2
  mean_expU = distributions[[distrib_idx]]$mean_expU
  
  distortion_functions = list(
    list(name = "Advanced_Linear",
         phi_X = function(U) 0.5 * U + log(2 / (exp(1) - exp(-1))),
         phi_Y = function(U) -0.5 * U - log(2 / (exp(1) - exp(-1)))),
    list(name = "Advanced_Quadratic",
         phi_X = function(U) U^2 - mean_U2 + 0.1 * exp(U - mean_expU),
         phi_Y = function(U) -(U^2 - mean_U2) - 0.1 * exp(U - mean_expU)),
    list(name = "Advanced_Periodic",
         phi_X = function(U) (1 + 0.5 * U) * sin(pi * U),
         phi_Y = function(U) -(1 + 0.5 * U) * sin(pi * U))
  )
  
  distrib_results = list()
  
  for (dist_idx in seq_along(distortion_functions)) {
    distortion_name = distortion_functions[[dist_idx]]$name
    phi_X_func = distortion_functions[[dist_idx]]$phi_X
    phi_Y_func = distortion_functions[[dist_idx]]$phi_Y
    
    results = matrix(NA, nrow = length(rho), ncol = 4)
    colnames(results) = c("Lower", "Upper", "Length", "Coverage")
    rownames(results) = as.character(rho)  # <--- Added this line
    
    
    for (ii in seq_along(rho)) {
      coverage = length = upper = lower = numeric(iter)
      
      for (jj in 1:iter) {
        sigma <- matrix(c(1, rho[ii], rho[ii], 1), 2)
        data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
        colnames(data) <- c("X", "Y")
        
        U = gen_U(n)
        phi_X = phi_X_func(U)
        phi_Y = phi_Y_func(U)
        xy.obs <- data + cbind(phi_X, phi_Y)
        
        e.X.est = e.Y.est = numeric(n)
        for (i in 1:n) {
          e.X.est[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
          e.Y.est[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
        }
        
        cov.e.est = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
        sig.e.X = mean(e.X.est^2) - mean(e.X.est)^2
        sig.e.Y = mean(e.Y.est^2) - mean(e.Y.est)^2
        rho.e.est = cov.e.est / sqrt(sig.e.X * sig.e.Y)
        
        rho.j = numeric(n)
        for (j in 1:n) {
          xy.jack = xy.obs[-j, ]
          U.jack = U[-j]
          
          Eix.j = Eiy.j = numeric(n - 1)
          for (k in 1:(n - 1)) {
            Eix.j[k] = xy.jack[k, 1] - log(get.NWK(xy.jack[, 1], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 1])))
            Eiy.j[k] = xy.jack[k, 2] - log(get.NWK(xy.jack[, 2], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 2])))
          }
          
          cov.j = mean(Eix.j * Eiy.j) - mean(Eix.j) * mean(Eiy.j)
          sig.jx = mean(Eix.j^2) - mean(Eix.j)^2
          sig.jy = mean(Eiy.j^2) - mean(Eiy.j)^2
          rho.jack = cov.j / sqrt(sig.jx * sig.jy)
          
          rho.j[j] = n * rho.e.est - (n - 1) * rho.jack
        }
        
        # AJEL step: add synthetic observation
        rho.j.adj = rho.j
        rho.j.adj[length(rho.j.adj) + 1] = -0.5 * log(n) * mean(rho.j.adj)
        
        il.ajel = el.test(rho.j.adj - rho[ii], 0)$'-2LLR'
        ci = findci(rho.j)
        
        upper[jj] = ci$Up
        lower[jj] = ci$Low
        length[jj] = ci$Up - ci$Low
        coverage[jj] = il.ajel < qchisq(0.975, 1)
      }
      
      results[ii, ] = c(mean(lower), mean(upper), mean(length), mean(coverage))
    }
    
    distrib_results[[distortion_name]] = results
  }
  
  all_results[[distrib_name]] = distrib_results
}

# Output
all_results


# MJEL
library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 50
rho = c(-0.9, -0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7, 0.9)
iter = 50

distributions = list(
  list(name = "Uniform", gen_U = function(n) runif(n, -1, 1), mean_U = 0, mean_U2 = 1/3, mean_expU = (exp(1) - exp(-1))/2),
  list(name = "Normal", gen_U = function(n) rnorm(n, 0, 1), mean_U = 0, mean_U2 = 1, mean_expU = exp(0.5)),
  list(name = "Beta", gen_U = function(n) rbeta(n, 2, 2) - 0.5, mean_U = 0, mean_U2 = 0.05, mean_expU = (exp(0.5) - exp(-0.5))/2)
)

all_results = list()

get.NWK <- function(x, u, small.u) {
  bw = density(u, kernel = c("epanechnikov"))$bw
  KX = K = numeric(length(u))
  for (j in 1:length(u)) {
    d = (u[j] - small.u) / bw
    weight = bw^(-1) * max(0, 0.75 * (1 - d^2))
    KX[j] = exp(x[j]) * weight
    K[j] = weight
  }
  return(sum(KX) / sum(K))
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    x.adj = x.vector - ci.val
    x.adj = x.adj - mean(x.adj)  # MJEL center
    el.test(x.adj, 0)
  }
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up = U, Low = L))
}

for (distrib_idx in seq_along(distributions)) {
  distrib_name = distributions[[distrib_idx]]$name
  gen_U = distributions[[distrib_idx]]$gen_U
  mean_U2 = distributions[[distrib_idx]]$mean_U2
  mean_expU = distributions[[distrib_idx]]$mean_expU
  
  distortion_functions = list(
    list(name = "Advanced_Linear",
         phi_X = function(U) 0.5 * U + log(2 / (exp(1) - exp(-1))),
         phi_Y = function(U) -0.5 * U - log(2 / (exp(1) - exp(-1)))),
    list(name = "Advanced_Quadratic",
         phi_X = function(U) U^2 - mean_U2 + 0.1 * exp(U - mean_expU),
         phi_Y = function(U) -(U^2 - mean_U2) - 0.1 * exp(U - mean_expU)),
    list(name = "Advanced_Periodic",
         phi_X = function(U) (1 + 0.5 * U) * sin(pi * U),
         phi_Y = function(U) -(1 + 0.5 * U) * sin(pi * U))
  )
  
  distrib_results = list()
  
  for (dist_idx in seq_along(distortion_functions)) {
    distortion_name = distortion_functions[[dist_idx]]$name
    phi_X_func = distortion_functions[[dist_idx]]$phi_X
    phi_Y_func = distortion_functions[[dist_idx]]$phi_Y
    
    results = matrix(NA, nrow = length(rho), ncol = 4)
    colnames(results) = c("Lower", "Upper", "Length", "Coverage")
    rownames(results) = as.character(rho)  # <--- Added this line
    
    
    for (ii in seq_along(rho)) {
      coverage = length = upper = lower = numeric(iter)
      
      for (jj in 1:iter) {
        sigma <- matrix(c(1, rho[ii], rho[ii], 1), 2)
        data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
        colnames(data) <- c("X", "Y")
        
        U = gen_U(n)
        phi_X = phi_X_func(U)
        phi_Y = phi_Y_func(U)
        xy.obs <- data + cbind(phi_X, phi_Y)
        
        e.X.est = e.Y.est = numeric(n)
        for (i in 1:n) {
          e.X.est[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
          e.Y.est[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
        }
        
        cov.e.est = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
        sig.e.X = mean(e.X.est^2) - mean(e.X.est)^2
        sig.e.Y = mean(e.Y.est^2) - mean(e.Y.est)^2
        rho.e.est = cov.e.est / sqrt(sig.e.X * sig.e.Y)
        
        rho.j = numeric(n)
        for (j in 1:n) {
          xy.jack = xy.obs[-j, ]
          U.jack = U[-j]
          
          Eix.j = Eiy.j = numeric(n - 1)
          for (k in 1:(n - 1)) {
            Eix.j[k] = xy.jack[k, 1] - log(get.NWK(xy.jack[, 1], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 1])))
            Eiy.j[k] = xy.jack[k, 2] - log(get.NWK(xy.jack[, 2], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 2])))
          }
          
          cov.j = mean(Eix.j * Eiy.j) - mean(Eix.j) * mean(Eiy.j)
          sig.jx = mean(Eix.j^2) - mean(Eix.j)^2
          sig.jy = mean(Eiy.j^2) - mean(Eiy.j)^2
          rho.jack = cov.j / sqrt(sig.jx * sig.jy)
          
          rho.j[j] = n * rho.e.est - (n - 1) * rho.jack
        }
        
        rho.j.centered = rho.j - mean(rho.j)
        
        il.mjel = el.test(rho.j.centered, 0)$'-2LLR'
        ci = findci(rho.j)
        
        upper[jj] = ci$Up
        lower[jj] = ci$Low
        length[jj] = ci$Up - ci$Low
        coverage[jj] = il.mjel < qchisq(0.975, 1)
      }
      
      results[ii, ] = c(mean(lower), mean(upper), mean(length), mean(coverage))
    }
    
    distrib_results[[distortion_name]] = results
  }
  
  all_results[[distrib_name]] = distrib_results
}

# Final output
all_results


# MAJEL
library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 50
rho = c(-0.9, -0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7, 0.9)
iter = 50

distributions = list(
  list(name = "Uniform", gen_U = function(n) runif(n, -1, 1), mean_U = 0, mean_U2 = 1/3, mean_expU = (exp(1) - exp(-1))/2),
  list(name = "Normal", gen_U = function(n) rnorm(n, 0, 1), mean_U = 0, mean_U2 = 1, mean_expU = exp(0.5)),
  list(name = "Beta", gen_U = function(n) rbeta(n, 2, 2) - 0.5, mean_U = 0, mean_U2 = 0.05, mean_expU = (exp(0.5) - exp(-0.5))/2)
)

all_results = list()

get.NWK <- function(x, u, small.u) {
  bw = density(u, kernel = c("epanechnikov"))$bw
  KX = K = numeric(length(u))
  for (j in 1:length(u)) {
    d = (u[j] - small.u) / bw
    weight = bw^(-1) * max(0, 0.75 * (1 - d^2))
    KX[j] = exp(x[j]) * weight
    K[j] = weight
  }
  return(sum(KX) / sum(K))
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    x.adj = x.vector - ci.val
    x.adj = x.adj - mean(x.adj)  # MAJEL center
    x.adj[length(x.adj) + 1] = -0.5 * log(length(x.adj)) * mean(x.adj)
    el.test(x.adj, 0)
  }
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up = U, Low = L))
}

for (distrib_idx in seq_along(distributions)) {
  distrib_name = distributions[[distrib_idx]]$name
  gen_U = distributions[[distrib_idx]]$gen_U
  mean_U2 = distributions[[distrib_idx]]$mean_U2
  mean_expU = distributions[[distrib_idx]]$mean_expU
  
  distortion_functions = list(
    list(name = "Advanced_Linear",
         phi_X = function(U) 0.5 * U + log(2 / (exp(1) - exp(-1))),
         phi_Y = function(U) -0.5 * U - log(2 / (exp(1) - exp(-1)))),
    list(name = "Advanced_Quadratic",
         phi_X = function(U) U^2 - mean_U2 + 0.1 * exp(U - mean_expU),
         phi_Y = function(U) -(U^2 - mean_U2) - 0.1 * exp(U - mean_expU)),
    list(name = "Advanced_Periodic",
         phi_X = function(U) (1 + 0.5 * U) * sin(pi * U),
         phi_Y = function(U) -(1 + 0.5 * U) * sin(pi * U))
  )
  
  distrib_results = list()
  
  for (dist_idx in seq_along(distortion_functions)) {
    distortion_name = distortion_functions[[dist_idx]]$name
    phi_X_func = distortion_functions[[dist_idx]]$phi_X
    phi_Y_func = distortion_functions[[dist_idx]]$phi_Y
    
    results = matrix(NA, nrow = length(rho), ncol = 4)
    colnames(results) = c("Lower", "Upper", "Length", "Coverage")
    rownames(results) = as.character(rho)  # <--- Added this line
    
    
    for (ii in seq_along(rho)) {
      coverage = length = upper = lower = numeric(iter)
      
      for (jj in 1:iter) {
        sigma <- matrix(c(1, rho[ii], rho[ii], 1), 2)
        data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
        colnames(data) <- c("X", "Y")
        
        U = gen_U(n)
        phi_X = phi_X_func(U)
        phi_Y = phi_Y_func(U)
        xy.obs <- data + cbind(phi_X, phi_Y)
        
        e.X.est = e.Y.est = numeric(n)
        for (i in 1:n) {
          e.X.est[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
          e.Y.est[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
        }
        
        cov.e.est = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
        sig.e.X = mean(e.X.est^2) - mean(e.X.est)^2
        sig.e.Y = mean(e.Y.est^2) - mean(e.Y.est)^2
        rho.e.est = cov.e.est / sqrt(sig.e.X * sig.e.Y)
        
        rho.j = numeric(n)
        for (j in 1:n) {
          xy.jack = xy.obs[-j, ]
          U.jack = U[-j]
          
          Eix.j = Eiy.j = numeric(n - 1)
          for (k in 1:(n - 1)) {
            Eix.j[k] = xy.jack[k, 1] - log(get.NWK(xy.jack[, 1], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 1])))
            Eiy.j[k] = xy.jack[k, 2] - log(get.NWK(xy.jack[, 2], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 2])))
          }
          
          cov.j = mean(Eix.j * Eiy.j) - mean(Eix.j) * mean(Eiy.j)
          sig.jx = mean(Eix.j^2) - mean(Eix.j)^2
          sig.jy = mean(Eiy.j^2) - mean(Eiy.j)^2
          rho.jack = cov.j / sqrt(sig.jx * sig.jy)
          
          rho.j[j] = n * rho.e.est - (n - 1) * rho.jack
        }
        
        rho.j.adj = rho.j - mean(rho.j)
        rho.j.adj[length(rho.j.adj) + 1] = -0.5 * log(n) * mean(rho.j.adj)
        
        il.majel = el.test(rho.j.adj, 0)$'-2LLR'
        ci = findci(rho.j)
        
        upper[jj] = ci$Up
        lower[jj] = ci$Low
        length[jj] = ci$Up - ci$Low
        coverage[jj] = il.majel < qchisq(0.975, 1)
      }
      
      results[ii, ] = c(mean(lower), mean(upper), mean(length), mean(coverage))
    }
    
    distrib_results[[distortion_name]] = results
  }
  
  all_results[[distrib_name]] = distrib_results
}

# Final output
all_results



# EL
library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 50
rho = c(-0.9, -0.5, 0, 0.5, 0.9)
iter = 50

distributions = list(
  list(name = "Uniform", gen_U = function(n) runif(n, -1, 1), mean_U = 0, mean_U2 = 1/3, mean_expU = (exp(1) - exp(-1))/2),
  list(name = "Normal", gen_U = function(n) rnorm(n, 0, 1), mean_U = 0, mean_U2 = 1, mean_expU = exp(0.5)),
  list(name = "Beta", gen_U = function(n) rbeta(n, 2, 2) - 0.5, mean_U = 0, mean_U2 = 0.05, mean_expU = (exp(0.5) - exp(-0.5))/2)
)

all_results = list()

get.NWK <- function(x, u, small.u) {
  bw = density(u, kernel = c("epanechnikov"))$bw
  KX = K = numeric(length(u))
  for (j in 1:length(u)) {
    d = (u[j] - small.u) / bw
    weight = bw^(-1) * max(0, 0.75 * (1 - d^2))
    KX[j] = exp(x[j]) * weight
    K[j] = weight
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

for (distrib_idx in seq_along(distributions)) {
  distrib_name = distributions[[distrib_idx]]$name
  gen_U = distributions[[distrib_idx]]$gen_U
  mean_U2 = distributions[[distrib_idx]]$mean_U2
  mean_expU = distributions[[distrib_idx]]$mean_expU
  
  distortion_functions = list(
    list(name = "Advanced_Linear",
         phi_X = function(U) 0.5 * U + log(2 / (exp(1) - exp(-1))),
         phi_Y = function(U) -0.5 * U - log(2 / (exp(1) - exp(-1)))),
    list(name = "Advanced_Quadratic",
         phi_X = function(U) U^2 - mean_U2 + 0.1 * exp(U - mean_expU),
         phi_Y = function(U) -(U^2 - mean_U2) - 0.1 * exp(U - mean_expU)),
    list(name = "Advanced_Periodic",
         phi_X = function(U) (1 + 0.5 * U) * sin(pi * U),
         phi_Y = function(U) -(1 + 0.5 * U) * sin(pi * U))
  )
  
  distrib_results = list()
  
  for (dist_idx in seq_along(distortion_functions)) {
    distortion_name = distortion_functions[[dist_idx]]$name
    phi_X_func = distortion_functions[[dist_idx]]$phi_X
    phi_Y_func = distortion_functions[[dist_idx]]$phi_Y
    
    results = matrix(NA, nrow = length(rho), ncol = 4)
    colnames(results) = c("Lower", "Upper", "Length", "Coverage")
    rownames(results) = as.character(rho)  # <--- Added this line
    
    
    for (ii in seq_along(rho)) {
      coverage = length = upper = lower = numeric(iter)
      
      for (jj in 1:iter) {
        sigma <- matrix(c(1, rho[ii], rho[ii], 1), 2)
        data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
        colnames(data) <- c("X", "Y")
        
        U = gen_U(n)
        phi_X = phi_X_func(U)
        phi_Y = phi_Y_func(U)
        xy.obs <- data + cbind(phi_X, phi_Y)
        
        e.X.est = e.Y.est = numeric(n)
        for (i in 1:n) {
          e.X.est[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i])) + log(mean(exp(xy.obs[, 1])))
          e.Y.est[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i])) + log(mean(exp(xy.obs[, 2])))
        }
        
        cov.e.est = mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
        sig.e.X = mean(e.X.est^2) - mean(e.X.est)^2
        sig.e.Y = mean(e.Y.est^2) - mean(e.Y.est)^2
        rho.e.est = cov.e.est / sqrt(sig.e.X * sig.e.Y)
        
        rho.j = numeric(n)
        for (j in 1:n) {
          xy.jack = xy.obs[-j, ]
          U.jack = U[-j]
          
          Eix.j = Eiy.j = numeric(n - 1)
          for (k in 1:(n - 1)) {
            Eix.j[k] = xy.jack[k, 1] - log(get.NWK(xy.jack[, 1], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 1])))
            Eiy.j[k] = xy.jack[k, 2] - log(get.NWK(xy.jack[, 2], U.jack, U.jack[k])) + log(mean(exp(xy.jack[, 2])))
          }
          
          cov.j = mean(Eix.j * Eiy.j) - mean(Eix.j) * mean(Eiy.j)
          sig.jx = mean(Eix.j^2) - mean(Eix.j)^2
          sig.jy = mean(Eiy.j^2) - mean(Eiy.j)^2
          rho.jack = cov.j / sqrt(sig.jx * sig.jy)
          
          rho.j[j] = n * rho.e.est - (n - 1) * rho.jack
        }
        
        il.el = el.test(rho.j - rho[ii], 0)$'-2LLR'
        ci = findci(rho.j)
        
        upper[jj] = ci$Up
        lower[jj] = ci$Low
        length[jj] = ci$Up - ci$Low
        coverage[jj] = il.el < qchisq(0.975, 1)
      }
      
      results[ii, ] = c(mean(lower), mean(upper), mean(length), mean(coverage))
    }
    
    distrib_results[[distortion_name]] = results
  }
  
  all_results[[distrib_name]] = distrib_results
}

# Final output
all_results
