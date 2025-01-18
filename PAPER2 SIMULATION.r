# Load necessary libraries
library(MASS)   # For multivariate normal data generation
library(emplik) # For empirical likelihood

# Parameters
n_values <- c(300, 500, 1000, 2000)  # Sample sizes
rho_values <- c(-0.9, -0.5, 0, 0.75) # Correlation coefficients
iterations <- 1000                   # Number of iterations

# Initialize results storage
results <- list()

define_bandwidth <- function(U, n) {
  # Silverman’s rule of thumb for bandwidth selection
  bw <- sd(U) * n^(-1/3)
  return(bw)
}

get.NWK <- function(x, u, small.u, bandwidth) {
  # Kernel smoothing function
  KX <- vector()
  K <- vector()
  for (j in 1:length(u)) {
    d <- (u[j] - small.u) / bandwidth
    KX[j] <- x[j] * max(0, 0.75 * (1 - d^2)) / bandwidth
    K[j] <- max(0, 0.75 * (1 - d^2)) / bandwidth
  }
  result <- sum(KX) / sum(K)
  return(result)
}

findci <- function(x.vector) {
  # Confidence interval using empirical likelihood
  cifunc <- function(ci.val, x.vector) {
    el.test(x.vector - ci.val, 0)
  }
  vectorMean <- mean(x.vector)
  L <- findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U <- findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Lower = L, Upper = U))
}

# Simulation Loop
for (n in n_values) {
  for (rho in rho_values) {
    coverage <- vector()
    ci_lengths <- vector()
    mse_values <- vector()

    for (iter in 1:iterations) {
      # 1. Generate True Data (X, Y)
      sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
      data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
      colnames(data) <- c("X", "Y")

      # 2. Generate Confounding Variable U
      U <- runif(n, 0, 1)

      # 3. Apply Distortion Functions
      psi_U <- 1.25 - 3 * (U - 0.5)^2
      phi_U <- 1 + 0.5 * cos(2 * pi * U)

      X_tilde <- psi_U * data[, "X"]
      Y_tilde <- phi_U * data[, "Y"]

      # 4. Residual-Based Estimators
      bw <- define_bandwidth(U, n)
      e_X_est <- sapply(1:n, function(j) X_tilde[j] / get.NWK(X_tilde, U, U[j], bw))
      e_Y_est <- sapply(1:n, function(j) Y_tilde[j] / get.NWK(Y_tilde, U, U[j], bw))

      cov_est <- mean(e_X_est * e_Y_est) - mean(e_X_est) * mean(e_Y_est)
      var_X <- mean(e_X_est^2) - mean(e_X_est)^2
      var_Y <- mean(e_Y_est^2) - mean(e_Y_est)^2
      rho_est <- cov_est / sqrt(var_X * var_Y)

      # 5. Confidence Interval Calculation
      rho_jackknife <- numeric(n)
      for (j in 1:n) {
        jack_X <- X_tilde[-j]
        jack_Y <- Y_tilde[-j]
        jack_U <- U[-j]

        jack_e_X <- sapply(1:(n-1), function(k) jack_X[k] / get.NWK(jack_X, jack_U, jack_U[k], bw))
        jack_e_Y <- sapply(1:(n-1), function(k) jack_Y[k] / get.NWK(jack_Y, jack_U, jack_U[k], bw))

        cov_jack <- mean(jack_e_X * jack_e_Y) - mean(jack_e_X) * mean(jack_e_Y)
        var_jack_X <- mean(jack_e_X^2) - mean(jack_e_X)^2
        var_jack_Y <- mean(jack_e_Y^2) - mean(jack_e_Y)^2

        rho_jackknife[j] <- n * rho_est - (n - 1) * (cov_jack / sqrt(var_jack_X * var_jack_Y))
      }

      ci <- findci(rho_jackknife)

      # 6. Evaluate Metrics
      coverage[iter] <- (ci$Lower <= rho) & (rho <= ci$Upper)
      ci_lengths[iter] <- ci$Upper - ci$Lower
      mse_values[iter] <- (rho_est - rho)^2
    }

    # Aggregate Results for this (n, rho) combination
    results[[paste0("n=", n, "_rho=", rho)]] <- list(
      Coverage = mean(coverage),
      CI_Length = mean(ci_lengths),
      MSE = mean(mse_values)
    )
  }
}

# Convert Results to Data Frame for Analysis
results_df <- do.call(rbind, lapply(names(results), function(key) {
  result <- results[[key]]
  data.frame(
    Scenario = key,
    Coverage = round(result$Coverage, 3),
    CI_Length = round(result$CI_Length, 3),
    MSE = round(result$MSE, 3)
  )
}))

# Display Results
print(results_df)

# Optional: Save Results to CSV
# write.csv(results_df, "paper2_simulation_results.csv", row.names = FALSE)
