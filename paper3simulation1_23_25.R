# Required Libraries
library(MASS) # For generating multivariate normal data
library(ks)   # For kernel density estimation

# Parameters for simulation
set.seed(123)
cor_values <- c(-0.9, -0.5, 0, 0.5, 0.9)  # Correlation coefficients
sample_sizes <- c(25, 50, 75, 100)        # Sample sizes
n_sim <- 10                             # Number of simulations

# Kernel function (Epanechnikov)
epanechnikov_kernel <- function(u) {
  k <- 0.75 * (1 - u^2)
  k[u > 1 | u < -1] <- 0
  return(k)
}

# Bandwidth calculation
calculate_bandwidth <- function(u) {
  sigma_u <- sd(u)
  n <- length(u)
  return(sigma_u * n^(-1/3))
}

# Residual-based estimator
calculate_residuals <- function(x_tilde, u, h) {
  n <- length(x_tilde)
  residuals <- numeric(n)
  for (i in 1:n) {
    weights <- epanechnikov_kernel((u - u[i]) / h)
    weights <- weights / sum(weights)
    residuals[i] <- x_tilde[i] - sum(weights * x_tilde)
  }
  return(residuals)
}

# Confidence Interval Calculation
calculate_ci <- function(rho_hat, method, sample_size) {
  if (method == "EL") {
    # Empirical Likelihood (placeholder)
    return(c(lower = rho_hat - 0.1, upper = rho_hat + 0.1))
  } else if (method == "JEL") {
    # Jackknife Empirical Likelihood (placeholder)
    return(c(lower = rho_hat - 0.12, upper = rho_hat + 0.12))
  } else if (method == "AJEL") {
    # Adjusted JEL (placeholder)
    return(c(lower = rho_hat - 0.15, upper = rho_hat + 0.15))
  } else if (method == "MJEL") {
    # Mean JEL (placeholder)
    return(c(lower = rho_hat - 0.14, upper = rho_hat + 0.14))
  } else if (method == "MAJEL") {
    # Mean Adjusted JEL (placeholder)
    return(c(lower = rho_hat - 0.16, upper = rho_hat + 0.16))
  }
}

# Simulation function
simulate <- function(n, rho, method, distribution = "uniform") {
  # Generate (X, Y)
  mu <- c(2, 4)
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  xy <- mvrnorm(n, mu, sigma)
  x <- xy[, 1]
  y <- xy[, 2]
  
  # Generate U based on the specified distribution
  if (distribution == "uniform") {
    u <- runif(n, 1, 7)
    psi_u <- u - 4
    phi_u <- 4 - u
  } else if (distribution == "normal") {
    u <- rnorm(n, 2, 1)
    psi_u <- u - 2
    phi_u <- 2 - u
  } else if (distribution == "beta") {
    u <- rbeta(n, 2, 8)
    psi_u <- u - 0.2
    phi_u <- 0.2 - u
  } else if (distribution == "weibull") {
    u <- rweibull(n, shape = 1.2, scale = 1)
    psi_u <- u - 0.9407
    phi_u <- 0.9407 - u
  }
  
  # Observed variables
  x_tilde <- x + psi_u
  y_tilde <- y + phi_u
  
  # Calculate residuals
  h <- calculate_bandwidth(u)
  residual_x <- calculate_residuals(x_tilde, u, h)
  residual_y <- calculate_residuals(y_tilde, u, h)
  
  # Residual-based estimator of correlation
  cov_xy <- mean(residual_x * residual_y) - mean(residual_x) * mean(residual_y)
  var_x <- mean(residual_x^2) - mean(residual_x)^2
  var_y <- mean(residual_y^2) - mean(residual_y)^2
  rho_hat <- cov_xy / sqrt(var_x * var_y)
  
  # Calculate confidence intervals
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
      
      # Calculate coverage probability and average length
      coverage <- mean(rho >= ci_lower & rho <= ci_upper)
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

# Save results to CSV
# write.csv(results, "simulation_results.csv", row.names = FALSE)
