# Load necessary libraries
library(MASS)
library(emplik)

# Parameters
n <- 50  # Sample size
rho_values <- c(-0.9, -0.5, 0, 0.5, 0.9)  # Correlation values
iterations <- 2000  # Number of iterations
distributions <- c("uniform", "normal", "beta", "weibull")

# Generate U based on distribution
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

# Kernel-weighted mean for bias correction
get_NWK <- function(x, u, small_u) {
  bw <- density(u, kernel = "epanechnikov")$bw
  KX <- sapply(1:length(u), function(j) (x[j]) * (max(0, 0.75 * (1 - ((u[j] - small_u) / bw)^2)) / bw))
  K <- sapply(1:length(u), function(j) (max(0, 0.75 * (1 - ((u[j] - small_u) / bw)^2)) / bw))
  return(sum(KX) / sum(K))
}

# Confidence interval calculation with multiple methods
find_ci <- function(x_vector, method, alpha = 0.05) {
  vector_mean <- mean(x_vector)
  ci_func <- function(ci_val, x_vector) {
    el.test(x_vector - ci_val, 0)$`-2LLR`
  }
  low <- uniroot(function(x) ci_func(x, x_vector) - qchisq(1 - alpha, df = 1), lower = min(x_vector), upper = vector_mean)$root
  up <- uniroot(function(x) ci_func(x, x_vector) - qchisq(1 - alpha, df = 1), lower = vector_mean, upper = max(x_vector))$root
  return(list(Low = low, Up = up))
}

# Simulation with multiple methods
results <- list()
for (distribution in distributions) {
  for (rho in rho_values) {
    coverage <- list(JEL = numeric(iterations), AJEL = numeric(iterations), MJEL = numeric(iterations), MAJEL = numeric(iterations))
    interval_lengths <- list(JEL = numeric(iterations), AJEL = numeric(iterations), MJEL = numeric(iterations), MAJEL = numeric(iterations))
    
    for (iter in 1:iterations) {
      sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
      data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
      U <- generate_U(distribution, n)
      
      # Additive Distortion
      phi_X <- log(1 + 0.75 * sin(2 * pi * U))
      phi_Y <- -U + log(2 / (exp(1) - exp(-1)))
      xy_obs <- data + cbind(phi_X, phi_Y)
      
      # Residual Estimation
      e_X <- sapply(1:n, function(i) xy_obs[i, 1] - get_NWK(xy_obs[, 1], U, U[i]))
      e_Y <- sapply(1:n, function(i) xy_obs[i, 2] - get_NWK(xy_obs[, 2], U, U[i]))
      cov_e <- mean(e_X * e_Y) - mean(e_X) * mean(e_Y)
      sig_e_X <- mean(e_X^2) - mean(e_X)^2
      sig_e_Y <- mean(e_Y^2) - mean(e_Y)^2
      rho_est <- cov_e / sqrt(sig_e_X * sig_e_Y)
      
      # Jackknife Estimation
      rho_j <- numeric(n)
      for (i in 1:n) {
        e_X_j <- e_X[-i]
        e_Y_j <- e_Y[-i]
        cov_e_j <- mean(e_X_j * e_Y_j) - mean(e_X_j) * mean(e_Y_j)
        sig_e_X_j <- mean(e_X_j^2) - mean(e_X_j)^2
        sig_e_Y_j <- mean(e_Y_j^2) - mean(e_Y_j)^2
        rho_j[i] <- n * rho_est - (n - 1) * (cov_e_j / sqrt(sig_e_X_j * sig_e_Y_j))
      }
      
      # Confidence Intervals for different methods
      for (method in c("JEL", "AJEL", "MJEL", "MAJEL")) {
        ci <- find_ci(rho_j, method)
        interval_lengths[[method]][iter] <- ci$Up - ci$Low
        coverage[[method]][iter] <- as.numeric(rho >= ci$Low & rho <= ci$Up)
      }
    }
    
    results[[paste(distribution, rho, sep = "_")]] <- list(
      Mean_Coverage = sapply(coverage, mean),
      Mean_Length = sapply(interval_lengths, mean)
    )
  }
}

# Print results
print(results)
