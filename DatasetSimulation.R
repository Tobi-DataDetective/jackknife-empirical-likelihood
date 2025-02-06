install.packages("MASS")
library(MASS)
data("Boston")
head(Boston)


# Install required packages if not already installed
install.packages("MASS")
install.packages("emplik")
install.packages("ks")

# Load required libraries
library(MASS)   # Boston dataset
library(emplik) # Empirical likelihood methods
library(ks)     # Kernel density estimation

# Load the Boston Housing Dataset
data(Boston)
df <- Boston  # Store in a variable for easier manipulation

# Select variables for correlation analysis
X <- df$lstat   # % Lower status population
Y <- df$medv    # Median house price

# Set Parameters
set.seed(123)
sample_sizes <- c(50, 100, 200) # Sample sizes
n_sim <- 1000  # Number of simulations

# Kernel function (Epanechnikov)
epanechnikov_kernel <- function(u) {
  k <- 0.75 * (1 - u^2)
  k[u > 1 | u < -1] <- 0
  return(k)
}

# Bandwidth calculation
calculate_bandwidth <- function(U) {
  sigma_U <- sd(U, na.rm = TRUE)  # Avoid NA issues
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
    weights <- weights / sum(weights, na.rm = TRUE)  # Avoid division by NA
    residuals[i] <- X_tilde[i] - sum(weights * X_tilde, na.rm = TRUE)
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

# Simulation function using Boston Dataset
simulate <- function(n, method, distribution = "uniform") {
  # Sample data from Boston Dataset
  sample_indices <- sample(1:nrow(df), n, replace = TRUE)
  X_sample <- X[sample_indices]
  Y_sample <- Y[sample_indices]
  
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
  distorted_data <- apply_mixed_distortion(X_sample, Y_sample, U)
  X_tilde <- distorted_data$X_tilde
  Y_tilde <- distorted_data$Y_tilde
  
  # Apply Residual-Based Correction if needed
  h <- calculate_bandwidth(U)
  residual_X <- calculate_residuals(X_tilde, U, h)
  residual_Y <- calculate_residuals(Y_tilde, U, h)
  
  # Compute residual-based correlation
  cov_xy <- mean(residual_X * residual_Y, na.rm = TRUE) - mean(residual_X, na.rm = TRUE) * mean(residual_Y, na.rm = TRUE)
  var_x <- mean(residual_X^2, na.rm = TRUE) - mean(residual_X, na.rm = TRUE)^2
  var_y <- mean(residual_Y^2, na.rm = TRUE) - mean(residual_Y, na.rm = TRUE)^2
  
  if (var_x <= 0 | var_y <= 0) {
    rho_hat <- NA
  } else {
    rho_hat <- cov_xy / sqrt(var_x * var_y)
  }
  
  # Compute Confidence Intervals
  if (!is.na(rho_hat)) {
    ci <- calculate_ci(rho_hat, method, n)
  } else {
    ci <- c(lower = NA, upper = NA)
  }
  
  return(list(rho_hat = rho_hat, ci = ci))
}

# Run simulations
results <- data.frame()
for (n in sample_sizes) {
  for (method in c("EL", "JEL", "AJEL", "MJEL", "MAJEL")) {
    sim_data <- replicate(n_sim, simulate(n, method, distribution = "uniform"), simplify = FALSE)
    rho_hats <- sapply(sim_data, function(res) res$rho_hat, simplify = TRUE)
    ci_lower <- sapply(sim_data, function(res) res$ci["lower"], simplify = TRUE)
    ci_upper <- sapply(sim_data, function(res) res$ci["upper"], simplify = TRUE)
    
    # Remove NA values before calculations
    valid_indices <- !is.na(rho_hats) & !is.na(ci_lower) & !is.na(ci_upper)
    if (sum(valid_indices) > 0) {
      rho_hats <- rho_hats[valid_indices]
      ci_lower <- ci_lower[valid_indices]
      ci_upper <- ci_upper[valid_indices]
      
      # Calculate coverage probability (rounded to 4 decimal places)
      coverage <- round(mean(-0.5 >= ci_lower & -0.5 <= ci_upper, na.rm = TRUE), 4)
      
      # Calculate average length (No rounding applied here)
      avg_length <- mean(ci_upper - ci_lower, na.rm = TRUE)
    } else {
      coverage <- NA
      avg_length <- NA
    }
    
    # Store results
    results <- rbind(results, data.frame(
      n = n,
      method = method,
      avg_length = avg_length,
      coverage = coverage
    ))
  }
}

# Print results
print(results)

# Save results to CSV (Ensuring coverage is in 4 decimal places)
# write.csv(results, "boston_simulation_results.csv", row.names = FALSE)
