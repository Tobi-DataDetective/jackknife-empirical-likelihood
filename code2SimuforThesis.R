# Install necessary packages (if not installed)
install.packages('emplik')
install.packages('kedd')

# Load libraries
library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS ########################
n = 50  # Sample size
rho = c(-0.9, -0.5, 0, 0.5, 0.9)  # Correlation values
iter = 1000  # Number of iterations

######## Define Matrix of Results ###################
coverage = matrix(NA, nrow = 5, ncol = 4)
lengths = matrix(NA, nrow = 5, ncol = 4)

########## FUNCTIONS #########################

# Kernel-weighted mean estimator
get_NWK = function(x, u, small_u) {
  bw = density(u, kernel = "epanechnikov")$bw
  KX = sapply(1:length(u), function(j) x[j] * max(0, 0.75 * (1 - ((u[j] - small_u) / bw)^2)))
  K = sapply(1:length(u), function(j) max(0, 0.75 * (1 - ((u[j] - small_u) / bw)^2)))
  return(sum(KX) / sum(K))
}

# Confidence interval estimation using Jackknife Empirical Likelihood
find_ci <- function(x_vector, method = "JEL") {
  ci_func <- function(ci_val, x_vector) {
    el.test(x_vector - ci_val, 0)$`-2LLR`
  }
  
  vector_mean <- mean(x_vector)
  lower <- uniroot(function(x) ci_func(x, x_vector) - qchisq(0.95, df = 1), lower = min(x_vector), upper = vector_mean)$root
  upper <- uniroot(function(x) ci_func(x, x_vector) - qchisq(0.95, df = 1), lower = vector_mean, upper = max(x_vector))$root
  
  if (method == "AJEL") {
    lower <- lower - 0.05 * abs(lower)
    upper <- upper + 0.05 * abs(upper)
  } else if (method == "MJEL") {
    lower <- lower - 0.1 * abs(lower)
    upper <- upper + 0.1 * abs(upper)
  } else if (method == "MAJEL") {
    lower <- lower - 0.15 * abs(lower)
    upper <- upper + 0.15 * abs(upper)
  }
  
  return(list(Low = lower, Up = upper))
}

############ ITERATIONS #######################
for (ii in 1:5) {
  for (jj in 1:iter) {
    
    ####### DATA GENERATION #####################
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma) 
    colnames(data) <- c("X", "Y")
    
    U = runif(n, -1, 1)  # Generate confounding variable
    
    ######## ADDITIVE DISTORTION #################
    phi_X = log(1 + 0.75 * sin(2 * pi * U))
    phi_Y = -U + log(2 / (exp(1) - exp(-1)))
    
    ###### OBSERVED DATA #########################
    xy_obs <- data + cbind(phi_X, phi_Y)
    
    ####### BIAS CORRECTED ESTIMATORS ##############
    e_X = sapply(1:n, function(i) xy_obs[i, 1] - get_NWK(xy_obs[, 1], U, U[i]))
    e_Y = sapply(1:n, function(i) xy_obs[i, 2] - get_NWK(xy_obs[, 2], U, U[i]))
    
    cov_e = mean(e_X * e_Y) - mean(e_X) * mean(e_Y)
    sig_e_X = mean(e_X^2) - mean(e_X)^2
    sig_e_Y = mean(e_Y^2) - mean(e_Y)^2
    
    rho_e_est = cov_e / sqrt(sig_e_X * sig_e_Y)
    
    ############ JACKKNIFE #######################
    rho_j = numeric(n)
    for (j in 1:n) {
      e_X_j = e_X[-j]
      e_Y_j = e_Y[-j]
      
      cov_e_j = mean(e_X_j * e_Y_j) - mean(e_X_j) * mean(e_Y_j)
      sig_e_X_j = mean(e_X_j^2) - mean(e_X_j)^2
      sig_e_Y_j = mean(e_Y_j^2) - mean(e_Y_j)^2
      
      rho_e_est_jack = cov_e_j / sqrt(sig_e_X_j * sig_e_Y_j)
      rho_j[j] = n * rho_e_est - (n - 1) * rho_e_est_jack
    }
    
    ######### JACKKNIFE EMPIRICAL LIKELIHOOD #################
    methods = c("JEL", "AJEL", "MJEL", "MAJEL")
    
    for (m in 1:length(methods)) {
      method = methods[m]
      ci = find_ci(rho_j, method)
      lengths[ii, m] = ci$Up - ci$Low
      coverage[ii, m] = mean(rho[ii] >= ci$Low & rho[ii] <= ci$Up)
    }
  }
}

# Combine results into a data frame
results = data.frame(
  Method = rep(methods, each = 5),
  Rho = rep(rho, times = length(methods)),
  Lower = coverage[, 1],
  Upper = coverage[, 2],
  Avg_Length = lengths[, 1],
  Coverage = coverage[, 3]
)

# Print results
print(results)
