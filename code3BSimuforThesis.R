# Install necessary packages
install.packages('emplik')

# Load required libraries
library(MASS)
library(emplik)

######### PARAMETERS ########################
n <- 50 # Sample size
rho_values <- c(-0.9, -0.5, 0, 0.5, 0.9) # Correlation coefficients
iterations <- 2000 # Number of iterations
distributions <- c("uniform", "normal", "beta", "weibull")

######### STORAGE ###########################
results <- list()

######### FUNCTIONS ##########################
# Generate U based on specified distribution
generate_U <- function(distribution, n) {
  if (distribution == "uniform") {
    return(runif(n, -1, 1))
  } else if (distribution == "normal") {
    return(rnorm(n, mean = 0, sd = 1))
  } else if (distribution == "beta") {
    return(rbeta(n, shape1 = 2, shape2 = 5))
  } else if (distribution == "weibull") {
    return(rweibull(n, shape = 1.2, scale = 1))
  }
}

# Kernel-weighted mean estimator for additive distortions
get_NWK <- function(x, u, small_u) {
  bw <- density(u, kernel = "epanechnikov")$bw
  KX <- sapply(1:length(u), function(j) (x[j]) * (max(0, 0.75 * (1 - ((u[j] - small_u) / bw)^2)) / bw))
  K <- sapply(1:length(u), function(j) max(0, 0.75 * (1 - ((u[j] - small_u) / bw)^2)) / bw)
  return(sum(KX) / sum(K))
}

# Confidence interval calculation
find_ci <- function(x_vector, alpha = 0.05) {
  vector_mean <- mean(x_vector)
  ci_func <- function(ci_val, x_vector) {
    el.test(x_vector - ci_val, 0)$-2LLR
  }
  low <- uniroot(function(x) ci_func(x, x_vector) - qchisq(1 - alpha, df = 1), lower = min(x_vector), upper = vector_mean)$root
  up <- uniroot(function(x) ci_func(x, x_vector) - qchisq(1 - alpha, df = 1), lower = vector_mean, upper = max(x_vector))$root
  return(list(Low = low, Up = up))
}

######### SIMULATION #########################
for (distribution in distributions) {
  for (rho in rho_values) {
    coverage <- numeric(iterations)
    interval_lengths <- numeric(iterations)
    
    for (iter in 1:iterations) {
      # Generate data
      sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
      data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
      U <- generate_U(distribution, n)
      
      # Additive distortions
      phi_X <- 0.5 * sin(2 * pi * U)
      phi_Y <- 0.75 * U - 0.3
      xy_obs <- data + cbind(phi_X, phi_Y)
      
      # Residual-based estimator
      e_X <- sapply(1:n, function(i) xy_obs[i, 1] - get_NWK(xy_obs[, 1], U, U[i]))
      e_Y <- sapply(1:n, function(i) xy_obs[i, 2] - get_NWK(xy_obs[, 2], U, U[i]))
      cov_e <- mean(e_X * e_Y) - mean(e_X) * mean(e_Y)
      sig_e_X <- mean(e_X^2) - mean(e_X)^2
      sig_e_Y <- mean(e_Y^2) - mean(e_Y)^2
      rho_est <- cov_e / sqrt(sig_e_X * sig_e_Y)
      
      # Jackknife estimation
      rho_j <- numeric(n)
      for (i in 1:n) {
        e_X_j <- e_X[-i]
        e_Y_j <- e_Y[-i]
        cov_e_j <- mean(e_X_j * e_Y_j) - mean(e_X_j) * mean(e_Y_j)
        sig_e_X_j <- mean(e_X_j^2) - mean(e_X_j)^2
        sig_e_Y_j <- mean(e_Y_j^2) - mean(e_Y_j)^2
        rho_j[i] <- n * rho_est - (n - 1) * (cov_e_j / sqrt(sig_e_X_j * sig_e_Y_j))
      }
      
      # Confidence intervals and coverage
      ci <- find_ci(rho_j)
      interval_lengths[iter] <- ci$Up - ci$Low
      coverage[iter] <- as.numeric(rho >= ci$Low & rho <= ci$Up)
    }
    
    # Store results
    results[[paste(distribution, rho, sep = "_")]] <- list(
      Mean_Coverage = mean(coverage),
      Mean_Length = mean(interval_lengths)
    )
  }
}

# Print results
results
