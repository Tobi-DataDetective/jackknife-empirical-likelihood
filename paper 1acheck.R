# Load necessary libraries
library(MASS)   # For multivariate normal data generation
library(emplik) # For empirical likelihood
library(kedd)   # For kernel density estimation
library(utils)  # For progress bar

# Parameters
n <- c(25, 50)  # Sample sizes
rho <- c(-0.9, 0, 0.9)  # Correlation coefficients
iter <- 100  # Adjusted iterations for faster testing; increase to 5000 for full runs
results <- list()  # To store simulation results
debug_results <- list()  # To store detailed debugging outputs

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
findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    el.test(x.vector - ci.val, 0)
  }
  vectorMean <- mean(x.vector)
  L <- findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U <- findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up = U, Low = L))
}

# Simulation
for (n_val in n) {
  for (rho_val in rho) {
    lower.jel <- vector()
    upper.jel <- vector()
    coverage.jel <- vector()
    length.jel <- vector()
    
    pb <- txtProgressBar(min = 0, max = iter, style = 3)
    
    for (i in 1:iter) {
      # Generate Data
      sigma <- matrix(c(1, rho_val, rho_val, 1), nrow = 2)
      data <- mvrnorm(n_val, mu = c(0, 0), Sigma = sigma)
      colnames(data) <- c("X", "Y")
      
      # Generate U (Confounding Variable)
      U <- runif(n_val, 0, 1)
      
      # Multiplicative Distortion Functions
      psi_U <- 1.25 - 3 * (U - 0.5)^2
      phi_U <- 1 + 0.5 * cos(2 * pi * U)
      
      # Observed Data
      X_tilde <- psi_U * data[, "X"]
      Y_tilde <- phi_U * data[, "Y"]
      
      # Residual-Based Estimators
      e.X.est <- sapply(1:n_val, function(j) {
        X_tilde[j] / get.NWK(X_tilde, U, U[j])
      })
      e.Y.est <- sapply(1:n_val, function(j) {
        Y_tilde[j] / get.NWK(Y_tilde, U, U[j])
      })
      
      cov.e.est <- mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
      sig.e.X <- mean(e.X.est^2) - mean(e.X.est)^2
      sig.e.Y <- mean(e.Y.est^2) - mean(e.Y.est)^2
      rho.e.est <- cov.e.est / sqrt(sig.e.X * sig.e.Y)
      
      # Jackknife Resampling
      rho.j <- numeric(n_val)
      for (j in 1:n_val) {
        jack_X <- X_tilde[-j]
        jack_Y <- Y_tilde[-j]
        jack_U <- U[-j]
        
        jack.e.X.est <- sapply(1:(n_val - 1), function(k) {
          jack_X[k] / get.NWK(jack_X, jack_U, jack_U[k])
        })
        jack.e.Y.est <- sapply(1:(n_val - 1), function(k) {
          jack_Y[k] / get.NWK(jack_Y, jack_U, jack_U[k])
        })
        
        cov.jack <- mean(jack.e.X.est * jack.e.Y.est) - mean(jack.e.X.est) * mean(jack.e.Y.est)
        sig.jack.X <- mean(jack.e.X.est^2) - mean(jack.e.X.est)^2
        sig.jack.Y <- mean(jack.e.Y.est^2) - mean(jack.e.Y.est)^2
        rho.j[j] <- n_val * rho.e.est - (n_val - 1) * (cov.jack / sqrt(sig.jack.X * sig.jack.Y))
      }
      
      # JEL Confidence Interval
      ci.jel <- findci(rho.j)
      lower.jel[i] <- ci.jel$Low
      upper.jel[i] <- ci.jel$Up
      length.jel[i] <- ci.jel$Up - ci.jel$Low
      coverage.jel[i] <- (ci.jel$Low <= rho_val) & (rho_val <= ci.jel$Up)
      
      # Store Debug Results
      debug_results[[paste0("n=", n_val, "_rho=", rho_val, "_iter=", i)]] <- list(
        Sigma = sigma,
        Psi_U = head(psi_U),
        Phi_U = head(phi_U),
        X_tilde = head(X_tilde),
        Y_tilde = head(Y_tilde),
        Estimated_Rho = rho.e.est,
        CI_Lower = ci.jel$Low,
        CI_Upper = ci.jel$Up,
        CI_Length = ci.jel$Up - ci.jel$Low
      )
      
      # Update progress bar
      setTxtProgressBar(pb, i)
    }
    
    close(pb)
    
    # Store Aggregate Results
    results[[paste0("n=", n_val, "_rho=", rho_val)]] <- list(
      Lower = mean(lower.jel),
      Upper = mean(upper.jel),
      Length = mean(length.jel),
      Coverage = mean(coverage.jel)
    )
  }
}

# Convert Debug Results to Tabular Format
debug_df <- do.call(rbind, lapply(debug_results, function(x) {
  data.frame(
    Sigma = paste(diag(x$Sigma), collapse = ", "),
    Psi_U = paste(round(x$Psi_U, 3), collapse = ", "),
    Phi_U = paste(round(x$Phi_U, 3), collapse = ", "),
    Estimated_Rho = round(x$Estimated_Rho, 3),
    CI_Lower = round(x$CI_Lower, 3),
    CI_Upper = round(x$CI_Upper, 3),
    CI_Length = round(x$CI_Length, 3)
  )
}))

# Display Results
head(debug_df)


# Convert Aggregate Results to Tabular Format
aggregate_results <- do.call(rbind, lapply(names(results), function(key) {
  result <- results[[key]]
  data.frame(
    Scenario = key,
    Lower = round(result$Lower, 3),
    Upper = round(result$Upper, 3),
    CI_Length = round(result$Length, 3),
    Coverage = round(result$Coverage, 3)
  )
}))

# Display the Aggregate Results
print(aggregate_results)

# Optional: Save Aggregate Results as a CSV for further analysis
# write.csv(aggregate_results, "aggregate_results.csv", row.names = FALSE)

