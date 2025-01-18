# Load necessary libraries
library(MASS)
library(emplik)
library(kedd)

# Parameters
n <- 50                         # Sample size
rho <- c(-0.9, -0.5, 0, 0.5, 0.9) # Correlation coefficients
iter <- 2000                    # Number of iterations

# Define Results Storage
coverage.jel <- vector()
length.jel <- vector()
upper.jel <- vector()
lower.jel <- vector()
results <- matrix(NA, nrow = length(rho), ncol = 4)

# Kernel Smoothing Function
get.NWK <- function(x, u, small.u) {
  bw <- density(u, kernel = "epanechnikov")$bw
  KX <- vector()
  K <- vector()
  for (j in 1:length(u)) {
    d <- (u[j] - small.u) / bw
    KX[j] <- x[j] * (bw^(-1) * max(0, 0.75 * (1 - d^2)))
    K[j] <- (bw^(-1) * max(0, 0.75 * (1 - d^2)))
  }
  result <- sum(KX) / sum(K)
  return(result)
}

# Confidence Interval Function
findci <- function(x.vector, AJEL = FALSE) {
  cifunc <- function(ci.val, x.vector, AJEL = AJEL) {
    if (AJEL) {
      x.vector <- x.vector - ci.val
      x.vector[length(x.vector) + 1] <- -0.5 * log(length(x.vector)) * mean(x.vector)
      el.test(x.vector, 0)
    } else {
      el.test(x.vector - ci.val, 0)
    }
  }
  vectorMean <- mean(x.vector)
  L <- findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector, AJEL = AJEL)$Low
  U <- findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector, AJEL = AJEL)$Up
  return(list(Up = U, Low = L))
}

# Simulation Loop
for (ii in 1:length(rho)) {
  for (jj in 1:iter) {
    
    # Data Generation
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
    colnames(data) <- c("X", "Y")
    
    # Generate U (Confounding Variable)
    U <- runif(n, -1, 1)
    
    # Multiplicative Distortion Functions
    psi_U <- 1.25 - 3 * (U - 0.5)^2
    phi_U <- 1 + 0.5 * cos(2 * pi * U)
    
    # Observed Data
    xy.obs <- data.frame(X = psi_U * data[, "X"], Y = phi_U * data[, "Y"])
    
    # Residual-Based Estimators
    e.Y.est <- sapply(1:n, function(i) {
      xy.obs$Y[i] / get.NWK(xy.obs$Y, U, U[i])
    })
    e.X.est <- sapply(1:n, function(i) {
      xy.obs$X[i] / get.NWK(xy.obs$X, U, U[i])
    })
    
    cov.e.est <- mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
    sig.e.X <- mean(e.X.est^2) - mean(e.X.est)^2
    sig.e.Y <- mean(e.Y.est^2) - mean(e.Y.est)^2
    rho.e.est <- cov.e.est / sqrt(sig.e.X * sig.e.Y)
    
    # Jackknifing
    rho.j <- numeric(n)
    for (j in 1:n) {
      jack.xy <- xy.obs[-j, ]
      jack.U <- U[-j]
      
      jack.e.Y.est <- sapply(1:(n - 1), function(k) {
        jack.xy$Y[k] / get.NWK(jack.xy$Y, jack.U, jack.U[k])
      })
      jack.e.X.est <- sapply(1:(n - 1), function(k) {
        jack.xy$X[k] / get.NWK(jack.xy$X, jack.U, jack.U[k])
      })
      
      cov.e.est.jack <- mean(jack.e.X.est * jack.e.Y.est) - mean(jack.e.X.est) * mean(jack.e.Y.est)
      sig.e.X.jack <- mean(jack.e.X.est^2) - mean(jack.e.X.est)^2
      sig.e.Y.jack <- mean(jack.e.Y.est^2) - mean(jack.e.Y.est)^2
      rho.j[j] <- n * rho.e.est - (n - 1) * (cov.e.est.jack / sqrt(sig.e.X.jack * sig.e.Y.jack))
    }
    
    # JEL Confidence Interval
    ci.jel <- findci(rho.j)
    upper.jel[jj] <- ci.jel$Up
    lower.jel[jj] <- ci.jel$Low
    length.jel[jj] <- upper.jel[jj] - lower.jel[jj]
    coverage.jel[jj] <- (ci.jel$Low <= rho[ii]) & (rho[ii] <= ci.jel$Up)
  }
  
  # Initialize a matrix to store results
  results <- matrix(NA, nrow = length(rho), ncol = 4)
  
  # After simulations, store the results for each rho
  for (ii in 1:length(rho)) {
    results[ii, ] <- c(mean(lower.jel), mean(upper.jel), mean(length.jel), mean(coverage.jel))
  }
  
  # Assign column and row names to the results matrix
  colnames(results) <- c("Lower Bound", "Upper Bound", "CI Length", "Coverage Probability")
  rownames(results) <- paste("Rho =", rho)
  
  # Convert results matrix to a data frame for better tabular output
  results_df <- as.data.frame(results)
  
  # Print results in a tabular format
  print(results_df)
  
  # Optional: Save results as a CSV file
  # write.csv(results_df, "simulation_results.csv", row.names = TRUE)
  
