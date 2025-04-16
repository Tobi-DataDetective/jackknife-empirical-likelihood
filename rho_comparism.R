#install.packages("emplik")
#install.packages("kedd")

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

get.NWK <- function(x, u, small.u, log_adjust = FALSE) {
  bw = density(u, kernel = "epanechnikov")$bw
  d = (u - small.u) / bw
  weight = pmax(0, 0.75 * (1 - d^2)) / bw
  if (log_adjust) {
    sum(exp(x) * weight) / sum(weight)  # For log-based estimator: E[exp(X)|U]
  } else {
    sum(x * weight) / sum(weight)  # Standard estimator: E[X|U]
  }
}

findci <- function(x.vector, method = "JEL") {
  cifunc <- function(ci.val, x.vector, method) {
    if (method == "AJEL") {
      x.adj = x.vector - ci.val
      n = length(x.vector)
      a_n = max(1, log(n) / 2)
      x.adj[n + 1] = -a_n * mean(x.adj) / n
      el.test(x.adj, 0)
    } else {
      el.test(x.vector - ci.val, 0)
    }
  }
  vectorMean <- mean(x.vector)
  L = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector, method = method)$Low
  U = findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector, method = method)$Up
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
    
    # Results matrices for each method and estimator
    methods = c("EL", "JEL", "AJEL")  # Removed MJEL and MAJEL
    results = list()
    for (method in methods) {
      results[[paste0(method, "_Std")]] = matrix(NA, nrow = length(rho), ncol = 4)
      results[[paste0(method, "_Log")]] = matrix(NA, nrow = length(rho), ncol = 4)
      colnames(results[[paste0(method, "_Std")]]) = colnames(results[[paste0(method, "_Log")]]) = c("Lower", "Upper", "Length", "Coverage")
      rownames(results[[paste0(method, "_Std")]]) = rownames(results[[paste0(method, "_Log")]]) = as.character(rho)
    }
    
    for (ii in seq_along(rho)) {
      coverage = length = upper = lower = list()
      for (method in methods) {
        coverage[[paste0(method, "_Std")]] = coverage[[paste0(method, "_Log")]] = numeric(iter)
        length[[paste0(method, "_Std")]] = length[[paste0(method, "_Log")]] = numeric(iter)
        upper[[paste0(method, "_Std")]] = upper[[paste0(method, "_Log")]] = numeric(iter)
        lower[[paste0(method, "_Std")]] = lower[[paste0(method, "_Log")]] = numeric(iter)
      }
      
      for (jj in 1:iter) {
        sigma <- matrix(c(1, rho[ii], rho[ii], 1), 2)
        data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
        colnames(data) <- c("X", "Y")
        
        U = gen_U(n)
        phi_X = phi_X_func(U)
        phi_Y = phi_Y_func(U)
        xy.obs <- data + cbind(phi_X, phi_Y)
        
        # Compute residuals using both estimators
        e.X.std = e.Y.std = e.X.log = e.Y.log = numeric(n)
        for (i in 1:n) {
          e.X.std[i] = xy.obs[i, 1] - get.NWK(xy.obs[, 1], U, U[i], log_adjust = FALSE)
          e.Y.std[i] = xy.obs[i, 2] - get.NWK(xy.obs[, 2], U, U[i], log_adjust = FALSE)
          e.X.log[i] = xy.obs[i, 1] - log(get.NWK(xy.obs[, 1], U, U[i], log_adjust = TRUE)) + log(mean(exp(xy.obs[, 1])))
          e.Y.log[i] = xy.obs[i, 2] - log(get.NWK(xy.obs[, 2], U, U[i], log_adjust = TRUE)) + log(mean(exp(xy.obs[, 2])))
        }
        
        # Correlation estimates
        cov.e.std = mean(e.X.std * e.Y.std) - mean(e.X.std) * mean(e.Y.std)
        sig.e.X.std = mean(e.X.std^2) - mean(e.X.std)^2
        sig.e.Y.std = mean(e.Y.std^2) - mean(e.Y.std)^2
        rho.e.std = cov.e.std / sqrt(sig.e.X.std * sig.e.Y.std)
        
        cov.e.log = mean(e.X.log * e.Y.log) - mean(e.X.log) * mean(e.Y.log)
        sig.e.X.log = mean(e.X.log^2) - mean(e.X.log)^2
        sig.e.Y.log = mean(e.Y.log^2) - mean(e.Y.log)^2
        rho.e.log = cov.e.log / sqrt(sig.e.X.log * sig.e.Y.log)
        
        # Jackknife pseudo-values for both estimators
        rho.j.std = rho.j.log = numeric(n)
        for (j in 1:n) {
          xy.jack = xy.obs[-j, ]
          U.jack = U[-j]
          Eix.j.std = Eiy.j.std = Eix.j.log = Eiy.j.log = numeric(n - 1)
          for (k in 1:(n - 1)) {
            Eix.j.std[k] = xy.jack[k, 1] - get.NWK(xy.jack[, 1], U.jack, U.jack[k], log_adjust = FALSE)
            Eiy.j.std[k] = xy.jack[k, 2] - get.NWK(xy.jack[, 2], U.jack, U.jack[k], log_adjust = FALSE)
            Eix.j.log[k] = xy.jack[k, 1] - log(get.NWK(xy.jack[, 1], U.jack, U.jack[k], log_adjust = TRUE)) + log(mean(exp(xy.jack[, 1])))
            Eiy.j.log[k] = xy.jack[k, 2] - log(get.NWK(xy.jack[, 2], U.jack, U.jack[k], log_adjust = TRUE)) + log(mean(exp(xy.jack[, 2])))
          }
          cov.j.std = mean(Eix.j.std * Eiy.j.std) - mean(Eix.j.std) * mean(Eiy.j.std)
          sig.jx.std = mean(Eix.j.std^2) - mean(Eix.j.std)^2
          sig.jy.std = mean(Eiy.j.std^2) - mean(Eiy.j.std)^2
          rho.jack.std = cov.j.std / sqrt(sig.jx.std * sig.jy.std)
          rho.j.std[j] = n * rho.e.std - (n - 1) * rho.jack.std
          
          cov.j.log = mean(Eix.j.log * Eiy.j.log) - mean(Eix.j.log) * mean(Eiy.j.log)
          sig.jx.log = mean(Eix.j.log^2) - mean(Eix.j.log)^2
          sig.jy.log = mean(Eiy.j.log^2) - mean(Eiy.j.log)^2
          rho.jack.log = cov.j.log / sqrt(sig.jx.log * sig.jy.log)
          rho.j.log[j] = n * rho.e.log - (n - 1) * rho.jack.log
        }
        
        # Compute confidence intervals for each method
        for (method in methods) {
          # Standard estimator
          ci = findci(rho.j.std, method = method)
          il = el.test(rho.j.std - rho[ii], 0)$'-2LLR'
          upper[[paste0(method, "_Std")]][jj] = ci$Up
          lower[[paste0(method, "_Std")]][jj] = ci$Low
          length[[paste0(method, "_Std")]][jj] = ci$Up - ci$Low
          coverage[[paste0(method, "_Std")]][jj] = il < qchisq(0.975, 1)
          
          # Log-adjusted estimator
          ci = findci(rho.j.log, method = method)
          il = el.test(rho.j.log - rho[ii], 0)$'-2LLR'
          upper[[paste0(method, "_Log")]][jj] = ci$Up
          lower[[paste0(method, "_Log")]][jj] = ci$Low
          length[[paste0(method, "_Log")]][jj] = ci$Up - ci$Low
          coverage[[paste0(method, "_Log")]][jj] = il < qchisq(0.975, 1)
        }
      }
      
      for (method in methods) {
        results[[paste0(method, "_Std")]][ii, ] = c(mean(lower[[paste0(method, "_Std")]]), mean(upper[[paste0(method, "_Std")]]), 
                                                    mean(length[[paste0(method, "_Std")]]), mean(coverage[[paste0(method, "_Std")]]))
        results[[paste0(method, "_Log")]][ii, ] = c(mean(lower[[paste0(method, "_Log")]]), mean(upper[[paste0(method, "_Log")]]), 
                                                    mean(length[[paste0(method, "_Log")]]), mean(coverage[[paste0(method, "_Log")]]))
      }
    }
    
    distrib_results[[distortion_name]] = results
  }
  
  all_results[[distrib_name]] = distrib_results
}

# Output results
library(knitr)
for (distrib_name in names(all_results)) {
  for (distortion_name in names(all_results[[distrib_name]])) {
    cat(sprintf("\n%s - %s\n", distrib_name, distortion_name))
    for (method in methods) {
      cat(sprintf("\n%s (Standard Estimator):\n", method))
      print(kable(all_results[[distrib_name]][[distortion_name]][[paste0(method, "_Std")]], digits = 3))
      cat(sprintf("\n%s (Log-Adjusted Estimator):\n", method))
      print(kable(all_results[[distrib_name]][[distortion_name]][[paste0(method, "_Log")]], digits = 3))
    }
  }
}


# 
# 
# Uniform - Advanced_Linear
# 
# EL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.973| -0.840|  0.133|     0.94|
#   |-0.5 | -0.754| -0.223|  0.531|     0.98|
#   |0    | -0.397|  0.296|  0.693|     0.98|
#   |0.5  |  0.226|  0.757|  0.531|     1.00|
#   |0.9  |  0.839|  0.969|  0.130|     0.96|
#   
#   EL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.986| -0.782|  0.204|     0.96|
#   |-0.5 | -0.744| -0.232|  0.512|     0.98|
#   |0    | -0.385|  0.292|  0.677|     0.96|
#   |0.5  |  0.218|  0.748|  0.530|     1.00|
#   |0.9  |  0.835|  0.970|  0.135|     1.00|
#   
#   JEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.973| -0.840|  0.133|     0.94|
#   |-0.5 | -0.754| -0.223|  0.531|     0.98|
#   |0    | -0.397|  0.296|  0.693|     0.98|
#   |0.5  |  0.226|  0.757|  0.531|     1.00|
#   |0.9  |  0.839|  0.969|  0.130|     0.96|
#   
#   JEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.986| -0.782|  0.204|     0.96|
#   |-0.5 | -0.744| -0.232|  0.512|     0.98|
#   |0    | -0.385|  0.292|  0.677|     0.96|
#   |0.5  |  0.218|  0.748|  0.530|     1.00|
#   |0.9  |  0.835|  0.970|  0.135|     1.00|
#   
#   AJEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.973| -0.839|  0.133|     0.94|
#   |-0.5 | -0.754| -0.223|  0.531|     0.98|
#   |0    | -0.397|  0.297|  0.694|     0.98|
#   |0.5  |  0.226|  0.757|  0.531|     1.00|
#   |0.9  |  0.839|  0.969|  0.130|     0.96|
#   
#   AJEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.986| -0.782|  0.204|     0.96|
#   |-0.5 | -0.744| -0.232|  0.512|     0.98|
#   |0    | -0.386|  0.292|  0.678|     0.96|
#   |0.5  |  0.218|  0.748|  0.530|     1.00|
#   |0.9  |  0.835|  0.970|  0.135|     1.00|
#   
#   Uniform - Advanced_Quadratic
# 
# EL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.966| -0.832|  0.134|     0.92|
#   |-0.5 | -0.764| -0.249|  0.515|     0.94|
#   |0    | -0.349|  0.339|  0.688|     0.98|
#   |0.5  |  0.319|  0.800|  0.481|     0.88|
#   |0.9  |  0.829|  0.969|  0.140|     0.94|
#   
#   EL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.985| -0.773|  0.212|     0.96|
#   |-0.5 | -0.748| -0.243|  0.505|     0.94|
#   |0    | -0.354|  0.313|  0.667|     0.96|
#   |0.5  |  0.305|  0.800|  0.494|     0.92|
#   |0.9  |  0.818|  0.970|  0.152|     0.94|
#   
#   JEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.966| -0.832|  0.134|     0.92|
#   |-0.5 | -0.764| -0.249|  0.515|     0.94|
#   |0    | -0.349|  0.339|  0.688|     0.98|
#   |0.5  |  0.319|  0.800|  0.481|     0.88|
#   |0.9  |  0.829|  0.969|  0.140|     0.94|
#   
#   JEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.985| -0.773|  0.212|     0.96|
#   |-0.5 | -0.748| -0.243|  0.505|     0.94|
#   |0    | -0.354|  0.313|  0.667|     0.96|
#   |0.5  |  0.305|  0.800|  0.494|     0.92|
#   |0.9  |  0.818|  0.970|  0.152|     0.94|
#   
#   AJEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.966| -0.832|  0.134|     0.92|
#   |-0.5 | -0.764| -0.249|  0.515|     0.94|
#   |0    | -0.350|  0.339|  0.688|     0.98|
#   |0.5  |  0.319|  0.800|  0.481|     0.88|
#   |0.9  |  0.829|  0.969|  0.140|     0.94|
#   
#   AJEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.985| -0.773|  0.213|     0.96|
#   |-0.5 | -0.748| -0.243|  0.505|     0.94|
#   |0    | -0.354|  0.313|  0.668|     0.96|
#   |0.5  |  0.305|  0.800|  0.495|     0.92|
#   |0.9  |  0.818|  0.970|  0.152|     0.94|
#   
#   Uniform - Advanced_Periodic
# 
# EL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.967| -0.831|  0.136|     0.96|
#   |-0.5 | -0.763| -0.265|  0.498|     1.00|
#   |0    | -0.302|  0.361|  0.663|     0.96|
#   |0.5  |  0.237|  0.761|  0.525|     0.98|
#   |0.9  |  0.820|  0.975|  0.155|     0.94|
#   
#   EL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.976| -0.749|  0.226|     0.96|
#   |-0.5 | -0.739| -0.260|  0.479|     0.98|
#   |0    | -0.303|  0.358|  0.660|     0.94|
#   |0.5  |  0.221|  0.756|  0.536|     0.96|
#   |0.9  |  0.810|  0.982|  0.172|     0.94|
#   
#   JEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.967| -0.831|  0.136|     0.96|
#   |-0.5 | -0.763| -0.265|  0.498|     1.00|
#   |0    | -0.302|  0.361|  0.663|     0.96|
#   |0.5  |  0.237|  0.761|  0.525|     0.98|
#   |0.9  |  0.820|  0.975|  0.155|     0.94|
#   
#   JEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.976| -0.749|  0.226|     0.96|
#   |-0.5 | -0.739| -0.260|  0.479|     0.98|
#   |0    | -0.303|  0.358|  0.660|     0.94|
#   |0.5  |  0.221|  0.756|  0.536|     0.96|
#   |0.9  |  0.810|  0.982|  0.172|     0.94|
#   
#   AJEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.967| -0.831|  0.136|     0.96|
#   |-0.5 | -0.763| -0.265|  0.498|     1.00|
#   |0    | -0.302|  0.361|  0.663|     0.96|
#   |0.5  |  0.237|  0.762|  0.525|     0.98|
#   |0.9  |  0.820|  0.975|  0.155|     0.94|
#   
#   AJEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.976| -0.749|  0.226|     0.96|
#   |-0.5 | -0.739| -0.260|  0.480|     0.98|
#   |0    | -0.303|  0.358|  0.661|     0.94|
#   |0.5  |  0.220|  0.756|  0.536|     0.96|
#   |0.9  |  0.810|  0.982|  0.172|     0.94|
#   
#   Normal - Advanced_Linear
# 
# EL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.977| -0.827|  0.150|     0.96|
#   |-0.5 | -0.788| -0.262|  0.526|     0.92|
#   |0    | -0.352|  0.377|  0.729|     0.98|
#   |0.5  |  0.229|  0.779|  0.550|     0.94|
#   |0.9  |  0.820|  0.976|  0.156|     0.96|
#   
#   EL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.975| -0.741|  0.234|     0.88|
#   |-0.5 | -0.761| -0.253|  0.509|     0.92|
#   |0    | -0.325|  0.393|  0.718|     0.98|
#   |0.5  |  0.227|  0.783|  0.555|     0.94|
#   |0.9  |  0.816|  0.980|  0.163|     0.98|
#   
#   JEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.977| -0.827|  0.150|     0.96|
#   |-0.5 | -0.788| -0.262|  0.526|     0.92|
#   |0    | -0.352|  0.377|  0.729|     0.98|
#   |0.5  |  0.229|  0.779|  0.550|     0.94|
#   |0.9  |  0.820|  0.976|  0.156|     0.96|
#   
#   JEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.975| -0.741|  0.234|     0.88|
#   |-0.5 | -0.761| -0.253|  0.509|     0.92|
#   |0    | -0.325|  0.393|  0.718|     0.98|
#   |0.5  |  0.227|  0.783|  0.555|     0.94|
#   |0.9  |  0.816|  0.980|  0.163|     0.98|
#   
#   AJEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.977| -0.827|  0.150|     0.96|
#   |-0.5 | -0.788| -0.262|  0.526|     0.92|
#   |0    | -0.352|  0.377|  0.729|     0.98|
#   |0.5  |  0.229|  0.779|  0.550|     0.94|
#   |0.9  |  0.820|  0.976|  0.156|     0.96|
#   
#   AJEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.975| -0.741|  0.234|     0.88|
#   |-0.5 | -0.762| -0.253|  0.509|     0.92|
#   |0    | -0.325|  0.393|  0.718|     0.98|
#   |0.5  |  0.227|  0.783|  0.555|     0.94|
#   |0.9  |  0.816|  0.980|  0.163|     0.98|
#   
#   Normal - Advanced_Quadratic
# 
# EL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.969| -0.825|  0.145|     0.98|
#   |-0.5 | -0.757| -0.241|  0.516|     0.98|
#   |0    | -0.401|  0.307|  0.708|     0.98|
#   |0.5  |  0.150|  0.743|  0.593|     0.94|
#   |0.9  |  0.787|  0.978|  0.191|     0.98|
#   
#   EL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.987| -0.733|  0.254|     0.88|
#   |-0.5 | -0.730| -0.212|  0.518|     1.00|
#   |0    | -0.385|  0.306|  0.691|     0.98|
#   |0.5  |  0.109|  0.734|  0.625|     0.90|
#   |0.9  |  0.749|  0.990|  0.242|     0.98|
#   
#   JEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.969| -0.825|  0.145|     0.98|
#   |-0.5 | -0.757| -0.241|  0.516|     0.98|
#   |0    | -0.401|  0.307|  0.708|     0.98|
#   |0.5  |  0.150|  0.743|  0.593|     0.94|
#   |0.9  |  0.787|  0.978|  0.191|     0.98|
#   
#   JEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.987| -0.733|  0.254|     0.88|
#   |-0.5 | -0.730| -0.212|  0.518|     1.00|
#   |0    | -0.385|  0.306|  0.691|     0.98|
#   |0.5  |  0.109|  0.734|  0.625|     0.90|
#   |0.9  |  0.749|  0.990|  0.242|     0.98|
#   
#   AJEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.970| -0.825|  0.145|     0.98|
#   |-0.5 | -0.757| -0.241|  0.516|     0.98|
#   |0    | -0.401|  0.307|  0.708|     0.98|
#   |0.5  |  0.149|  0.743|  0.594|     0.94|
#   |0.9  |  0.787|  0.978|  0.191|     0.98|
#   
#   AJEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.987| -0.733|  0.254|     0.88|
#   |-0.5 | -0.730| -0.212|  0.518|     1.00|
#   |0    | -0.385|  0.307|  0.692|     0.98|
#   |0.5  |  0.108|  0.734|  0.626|     0.90|
#   |0.9  |  0.748|  0.991|  0.242|     0.98|
#   
#   Normal - Advanced_Periodic
# 
# EL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.968| -0.827|  0.141|     0.98|
#   |-0.5 | -0.764| -0.259|  0.505|     1.00|
#   |0    | -0.421|  0.305|  0.726|     0.98|
#   |0.5  |  0.192|  0.777|  0.586|     0.98|
#   |0.9  |  0.741|  0.989|  0.248|     0.94|
#   
#   EL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -1.007| -0.742|  0.265|     0.88|
#   |-0.5 | -0.733| -0.227|  0.505|     0.92|
#   |0    | -0.421|  0.308|  0.729|     0.96|
#   |0.5  |  0.145|  0.751|  0.606|     0.98|
#   |0.9  |  0.711|  1.017|  0.306|     0.92|
#   
#   JEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.968| -0.827|  0.141|     0.98|
#   |-0.5 | -0.764| -0.259|  0.505|     1.00|
#   |0    | -0.421|  0.305|  0.726|     0.98|
#   |0.5  |  0.192|  0.777|  0.586|     0.98|
#   |0.9  |  0.741|  0.989|  0.248|     0.94|
#   
#   JEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -1.007| -0.742|  0.265|     0.88|
#   |-0.5 | -0.733| -0.227|  0.505|     0.92|
#   |0    | -0.421|  0.308|  0.729|     0.96|
#   |0.5  |  0.145|  0.751|  0.606|     0.98|
#   |0.9  |  0.711|  1.017|  0.306|     0.92|
#   
#   AJEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.968| -0.827|  0.141|     0.98|
#   |-0.5 | -0.765| -0.259|  0.506|     1.00|
#   |0    | -0.421|  0.305|  0.726|     0.98|
#   |0.5  |  0.191|  0.777|  0.586|     0.98|
#   |0.9  |  0.741|  0.989|  0.249|     0.94|
#   
#   AJEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -1.007| -0.742|  0.266|     0.88|
#   |-0.5 | -0.733| -0.227|  0.506|     0.92|
#   |0    | -0.422|  0.308|  0.730|     0.96|
#   |0.5  |  0.144|  0.751|  0.607|     0.98|
#   |0.9  |  0.711|  1.017|  0.306|     0.92|
#   
#   Beta - Advanced_Linear
# 
# EL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.973| -0.835|  0.139|     0.90|
#   |-0.5 | -0.780| -0.278|  0.502|     0.94|
#   |0    | -0.329|  0.352|  0.681|     0.98|
#   |0.5  |  0.218|  0.745|  0.527|     0.96|
#   |0.9  |  0.831|  0.975|  0.144|     0.98|
#   
#   EL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -1.002| -0.788|  0.213|     0.96|
#   |-0.5 | -0.745| -0.271|  0.474|     0.94|
#   |0    | -0.321|  0.344|  0.665|     0.98|
#   |0.5  |  0.203|  0.741|  0.538|     0.98|
#   |0.9  |  0.827|  0.979|  0.152|     1.00|
#   
#   JEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.973| -0.835|  0.139|     0.90|
#   |-0.5 | -0.780| -0.278|  0.502|     0.94|
#   |0    | -0.329|  0.352|  0.681|     0.98|
#   |0.5  |  0.218|  0.745|  0.527|     0.96|
#   |0.9  |  0.831|  0.975|  0.144|     0.98|
#   
#   JEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -1.002| -0.788|  0.213|     0.96|
#   |-0.5 | -0.745| -0.271|  0.474|     0.94|
#   |0    | -0.321|  0.344|  0.665|     0.98|
#   |0.5  |  0.203|  0.741|  0.538|     0.98|
#   |0.9  |  0.827|  0.979|  0.152|     1.00|
#   
#   AJEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.974| -0.835|  0.139|     0.90|
#   |-0.5 | -0.780| -0.278|  0.502|     0.94|
#   |0    | -0.329|  0.353|  0.682|     0.98|
#   |0.5  |  0.217|  0.745|  0.528|     0.96|
#   |0.9  |  0.831|  0.975|  0.144|     0.98|
#   
#   AJEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -1.002| -0.788|  0.213|     0.96|
#   |-0.5 | -0.745| -0.270|  0.474|     0.94|
#   |0    | -0.321|  0.344|  0.665|     0.98|
#   |0.5  |  0.203|  0.741|  0.538|     0.98|
#   |0.9  |  0.827|  0.979|  0.152|     1.00|
#   
#   Beta - Advanced_Quadratic
# 
# EL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.977| -0.839|  0.138|     0.98|
#   |-0.5 | -0.743| -0.231|  0.512|     0.96|
#   |0    | -0.332|  0.334|  0.666|     0.96|
#   |0.5  |  0.220|  0.749|  0.528|     0.94|
#   |0.9  |  0.839|  0.974|  0.135|     0.96|
#   
#   EL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.977| -0.730|  0.247|     0.96|
#   |-0.5 | -0.717| -0.218|  0.498|     0.98|
#   |0    | -0.335|  0.327|  0.661|     0.96|
#   |0.5  |  0.206|  0.735|  0.529|     0.94|
#   |0.9  |  0.833|  0.978|  0.145|     0.94|
#   
#   JEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.977| -0.839|  0.138|     0.98|
#   |-0.5 | -0.743| -0.231|  0.512|     0.96|
#   |0    | -0.332|  0.334|  0.666|     0.96|
#   |0.5  |  0.220|  0.749|  0.528|     0.94|
#   |0.9  |  0.839|  0.974|  0.135|     0.96|
#   
#   JEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.977| -0.730|  0.247|     0.96|
#   |-0.5 | -0.717| -0.218|  0.498|     0.98|
#   |0    | -0.335|  0.327|  0.661|     0.96|
#   |0.5  |  0.206|  0.735|  0.529|     0.94|
#   |0.9  |  0.833|  0.978|  0.145|     0.94|
#   
#   AJEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.977| -0.839|  0.138|     0.98|
#   |-0.5 | -0.743| -0.231|  0.512|     0.96|
#   |0    | -0.332|  0.334|  0.666|     0.96|
#   |0.5  |  0.220|  0.749|  0.529|     0.94|
#   |0.9  |  0.839|  0.974|  0.135|     0.96|
#   
#   AJEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.977| -0.730|  0.247|     0.96|
#   |-0.5 | -0.717| -0.218|  0.499|     0.98|
#   |0    | -0.335|  0.327|  0.662|     0.96|
#   |0.5  |  0.206|  0.735|  0.530|     0.94|
#   |0.9  |  0.833|  0.978|  0.145|     0.94|
#   
#   Beta - Advanced_Periodic
# 
# EL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.966| -0.835|  0.131|     0.96|
#   |-0.5 | -0.811| -0.272|  0.538|     0.96|
#   |0    | -0.401|  0.302|  0.703|     0.92|
#   |0.5  |  0.227|  0.758|  0.531|     0.98|
#   |0.9  |  0.834|  0.970|  0.136|     0.96|
#   
#   EL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.977| -0.773|  0.205|     0.96|
#   |-0.5 | -0.780| -0.266|  0.514|     0.94|
#   |0    | -0.387|  0.290|  0.678|     0.90|
#   |0.5  |  0.208|  0.749|  0.540|     1.00|
#   |0.9  |  0.825|  0.970|  0.146|     0.96|
#   
#   JEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.966| -0.835|  0.131|     0.96|
#   |-0.5 | -0.811| -0.272|  0.538|     0.96|
#   |0    | -0.401|  0.302|  0.703|     0.92|
#   |0.5  |  0.227|  0.758|  0.531|     0.98|
#   |0.9  |  0.834|  0.970|  0.136|     0.96|
#   
#   JEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.977| -0.773|  0.205|     0.96|
#   |-0.5 | -0.780| -0.266|  0.514|     0.94|
#   |0    | -0.387|  0.290|  0.678|     0.90|
#   |0.5  |  0.208|  0.749|  0.540|     1.00|
#   |0.9  |  0.825|  0.970|  0.146|     0.96|
#   
#   AJEL (Standard Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.966| -0.835|  0.131|     0.96|
#   |-0.5 | -0.811| -0.272|  0.539|     0.96|
#   |0    | -0.401|  0.302|  0.703|     0.92|
#   |0.5  |  0.227|  0.758|  0.532|     0.98|
#   |0.9  |  0.834|  0.970|  0.136|     0.96|
#   
#   AJEL (Log-Adjusted Estimator):
#   
#   
#   |     |  Lower|  Upper| Length| Coverage|
#   |:----|------:|------:|------:|--------:|
#   |-0.9 | -0.978| -0.773|  0.205|     0.96|
#   |-0.5 | -0.781| -0.266|  0.515|     0.94|
#   |0    | -0.388|  0.291|  0.678|     0.90|
#   |0.5  |  0.208|  0.749|  0.541|     1.00|
#   |0.9  |  0.825|  0.971|  0.146|     0.96|