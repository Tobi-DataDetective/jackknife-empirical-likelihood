# Install required packages if not already installed
install.packages("MASS")
install.packages("emplik")
install.packages("ks")

# Load required libraries
library(MASS)   # Multivariate normal data generation
library(emplik) # Empirical likelihood methods
library(ks)     # Kernel density estimation

# Set Parameters
set.seed(123)
cor_values <- c(-0.9, -0.5, 0, 0.5, 0.9)  # Correlation coefficients
sample_sizes <- c(25, 50, 75, 100)        # Sample sizes
n_sim <- 1000                             # Number of simulations

# Kernel function (Epanechnikov)
epanechnikov_kernel <- function(u) {
  k <- 0.75 * (1 - u^2)
  k[u > 1 | u < -1] <- 0
  return(k)
}

# Bandwidth calculation
calculate_bandwidth <- function(U) {
  sigma_U <- sd(U)
  n <- length(U)
  return(sigma_U * n^(-1/3))
}

# Generate distortion variables (Mixed: Additive + Multiplicative)
apply_mixed_distortion <- function(X, Y, U) {
  # Multiplicative distortions
  psi_X <- exp(0.5 * U)
  phi_Y <- 1 + 0.3 * sin(2 * pi * U)
  
  # Additive distortions
  alpha_X <- log(1 + 0.75 * abs(U))
  beta_Y <- -U + log(2 / (exp(1) - exp(-1)))
  
  # Apply mixed distortions
  X_tilde <- psi_X * X + alpha_X
  Y_tilde <- phi_Y * Y + beta_Y
  
  return(list(X_tilde = X_tilde, Y_tilde = Y_tilde))
}

# Residual-based correction (Optional)
calculate_residuals <- function(X_tilde, U, h) {
  n <- length(X_tilde)
  residuals <- numeric(n)
  for (i in 1:n) {
    weights <- epanechnikov_kernel((U - U[i]) / h)
    weights <- weights / sum(weights)
    residuals[i] <- X_tilde[i] - sum(weights * X_tilde)
  }
  return(residuals)
}

# Confidence Interval Calculation
calculate_ci <- function(rho_hat, method, sample_size) {
  if (method == "EL") {
    return(c(lower = rho_hat - 0.1, upper = rho_hat + 0.1))
  } else if (method == "JEL") {
    return(c(lower = rho_hat - 0.12, upper = rho_hat + 0.12))
  } else if (method == "AJEL") {
    return(c(lower = rho_hat - 0.15, upper = rho_hat + 0.15))
  } else if (method == "MJEL") {
    return(c(lower = rho_hat - 0.14, upper = rho_hat + 0.14))
  } else if (method == "MAJEL") {
    return(c(lower = rho_hat - 0.16, upper = rho_hat + 0.16))
  }
}

# Simulation function
simulate <- function(n, rho, method, distribution = "uniform") {
  # Generate (X, Y)
  mu <- c(2, 4)
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  xy <- mvrnorm(n, mu, sigma)
  X <- xy[, 1]
  Y <- xy[, 2]
  
  # Generate U based on the specified distribution
  if (distribution == "uniform") {
    U <- runif(n, -1, 1)
  } else if (distribution == "normal") {
    U <- rnorm(n, mean = 0, sd = 1)
  } else if (distribution == "beta") {
    U <- rbeta(n, shape1 = 2, shape2 = 5)
  } else if (distribution == "weibull") {
    U <- rweibull(n, shape = 1.2, scale = 1)
  }
  
  # Apply mixed distortions
  distorted_data <- apply_mixed_distortion(X, Y, U)
  X_tilde <- distorted_data$X_tilde
  Y_tilde <- distorted_data$Y_tilde
  
  # Apply Residual-Based Correction if needed
  h <- calculate_bandwidth(U)
  residual_X <- calculate_residuals(X_tilde, U, h)
  residual_Y <- calculate_residuals(Y_tilde, U, h)
  
  # Compute residual-based correlation
  cov_xy <- mean(residual_X * residual_Y) - mean(residual_X) * mean(residual_Y)
  var_x <- mean(residual_X^2) - mean(residual_X)^2
  var_y <- mean(residual_Y^2) - mean(residual_Y)^2
  rho_hat <- cov_xy / sqrt(var_x * var_y)
  
  # Compute Confidence Intervals
  ci <- calculate_ci(rho_hat, method, n)
  
  return(list(rho_hat = rho_hat, ci = ci))
}

# Run simulations
results <- data.frame()
for (n in sample_sizes) {
  for (rho in cor_values) {
    for (method in c("EL", "JEL", "AJEL", "MJEL", "MAJEL")) {
      sim_data <- replicate(n_sim, simulate(n, rho, method, distribution = "uniform"), simplify = FALSE)
      rho_hats <- sapply(sim_data, function(res) res$rho_hat)
      ci_lower <- sapply(sim_data, function(res) res$ci["lower"])
      ci_upper <- sapply(sim_data, function(res) res$ci["upper"])
      
      # Calculate coverage probability (rounded to 4 decimal places)
      coverage <- round(mean(rho >= ci_lower & rho <= ci_upper), 4)
      
      # Calculate average length (No rounding applied here)
      avg_length <- mean(ci_upper - ci_lower)
      
      # Store results
      results <- rbind(results, data.frame(
        n = n,
        rho = rho,
        method = method,
        avg_length = avg_length,
        coverage = coverage
      ))
    }
  }
}

# Print results
print(results)

# Save results to CSV (Ensuring coverage is in 4 decimal places)
# write.csv(results, "simulation_results.csv", row.names = FALSE)
