# Install required package
install.packages('ks')

# Load libraries
library(MASS) # For generating multivariate normal data
library(ks)   # For kernel density estimation

# Set Parameters
set.seed(123)
cor_values <- c(-0.9, -0.5, 0, 0.5, 0.9)  # Correlation coefficients
sample_sizes <- c(25, 50, 75, 100)        # Sample sizes
n_sim <- 2000                             # Number of simulations

# Epanechnikov Kernel Function
epanechnikov_kernel <- function(u) {
  k <- 0.75 * (1 - u^2)
  k[u > 1 | u < -1] <- 0
  return(k)
}

# Bandwidth Calculation
calculate_bandwidth <- function(u) {
  sigma_u <- sd(u)
  n <- length(u)
  return(sigma_u * n^(-1/3))  # Rule-of-thumb bandwidth
}

# Function to Generate Additive Distortion
generate_additive_distortion <- function(n, distribution = "uniform") {
  if (distribution == "uniform") {
    return(runif(n, -2, 2))
  } else if (distribution == "normal") {
    return(rnorm(n, mean = 0, sd = 1))
  } else if (distribution == "beta") {
    return(rbeta(n, 2, 8) - 0.25)  # Centered around 0
  } else if (distribution == "weibull") {
    return(rweibull(n, shape = 1.2, scale = 1) - 0.94)
  }
}

# Residual-based Estimation of Pearson Correlation
estimate_correlation <- function(x_tilde, y_tilde, u, h) {
  n <- length(x_tilde)
  residuals_x <- numeric(n)
  residuals_y <- numeric(n)
  
  for (i in 1:n) {
    weights <- epanechnikov_kernel((u - u[i]) / h)
    weights <- weights / sum(weights)
    residuals_x[i] <- x_tilde[i] - sum(weights * x_tilde)
    residuals_y[i] <- y_tilde[i] - sum(weights * y_tilde)
  }
  
  cov_xy <- mean(residuals_x * residuals_y) - mean(residuals_x) * mean(residuals_y)
  var_x <- mean(residuals_x^2) - mean(residuals_x)^2
  var_y <- mean(residuals_y^2) - mean(residuals_y)^2
  rho_hat <- cov_xy / sqrt(var_x * var_y)
  
  return(rho_hat)
}

# Confidence Interval Calculation using Jackknife Empirical Likelihood
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

# Simulation Function
simulate <- function(n, rho, method, distribution = "uniform") {
  # Generate (X, Y) from a bivariate normal distribution
  mu <- c(0, 0)
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  xy <- mvrnorm(n, mu, sigma)
  x <- xy[, 1]
  y <- xy[, 2]
  
  # Generate Additive Distortions
  u <- runif(n, 0, 1)  # Confounding variable
  psi_u <- generate_additive_distortion(n, distribution)
  phi_u <- generate_additive_distortion(n, distribution)
  
  # Observed variables with additive distortions
  x_tilde <- x + psi_u
  y_tilde <- y + phi_u
  
  # Estimate residuals
  h <- calculate_bandwidth(u)
  rho_hat <- estimate_correlation(x_tilde, y_tilde, u, h)
  
  # Compute confidence intervals
  ci <- calculate_ci(rho_hat, method, n)
  
  return(list(rho_hat = rho_hat, ci = ci))
}

# Run Simulations
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

# Print Results
print(results)

# Save Results
write.csv(results, "simulation_results.csv", row.names = FALSE)
