library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################

n = 75
rho = c(-0.9, -0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7, 0.9)
iter = 100

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
