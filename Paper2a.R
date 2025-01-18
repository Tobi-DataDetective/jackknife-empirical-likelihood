library(MASS)
library(emplik)
library(kedd)

######### PARAMETERS########################

n_values <- c(50, 100, 200, 300)  # Reduced for testing
rho <- c(-0.9, -0.5, 0, 0.75)
iter <- 100  # Reduced iterations for testing

########Define Matrix of Results###################
results <- matrix(NA, nrow = length(rho), ncol = 4)
colnames(results) <- c("Lower Bound", "Upper Bound", "CI Length", "Coverage Probability")
rownames(results) <- paste("Rho =", rho)

########## FUNCTIONS#########################

get.NWK <- function(x, u, small.u) {
  bw <- density(u, kernel = c("epanechnikov"))$bw
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

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    el.test(x.vector - ci.val, 0)
  }
  vectorMean <- mean(x.vector)
  L <- findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U <- findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Up = U, Low = L))
}

############SIMULATION#######################

for (ii in seq_along(rho)) {
  
  lower.jel <- vector()
  upper.jel <- vector()
  coverage.jel <- vector()
  length.jel <- vector()
  
  for (jj in 1:iter) {
    
    ####### DATA GENERATION#####################
    
    sigma <- matrix(c(1, rho[ii], rho[ii], 1), nrow = 2)
    data <- mvrnorm(n_values[ii %% length(n_values)], mu = c(0, 0), Sigma = sigma)
    colnames(data) <- c("X", "Y")
    
    U <- runif(n_values[ii %% length(n_values)], 0, 1)  # Adjusted U to Uniform(0, 1)
    
    ########DISTORTING FUNCTION#################
    
    psi_U <- 1.25 - 3 * (U - 0.5)^2
    phi_U <- 1 + 0.5 * cos(2 * pi * U)
    
    ###### OBSERVED DATA#########################
    
    X_tilde <- psi_U * data[, "X"]
    Y_tilde <- phi_U * data[, "Y"]
    
    #######CALCULATING ESTIMATORS##############
    
    e.X.est <- sapply(1:length(U), function(j) X_tilde[j] / get.NWK(X_tilde, U, U[j]))
    e.Y.est <- sapply(1:length(U), function(j) Y_tilde[j] / get.NWK(Y_tilde, U, U[j]))
    
    cov.e.est <- mean(e.X.est * e.Y.est) - mean(e.X.est) * mean(e.Y.est)
    sig.e.X <- mean(e.X.est^2) - mean(e.X.est)^2
    sig.e.Y <- mean(e.Y.est^2) - mean(e.Y.est)^2
    rho.e.est <- cov.e.est / sqrt(sig.e.X * sig.e.Y)
    
    ############JACKKNIFING#######################
    
    rho.j <- numeric(length(U))
    for (j in 1:length(U)) {
      jack_X <- X_tilde[-j]
      jack_Y <- Y_tilde[-j]
      jack_U <- U[-j]
      
      jack.e.X.est <- sapply(1:(length(U) - 1), function(k) {
        jack_X[k] / get.NWK(jack_X, jack_U, jack_U[k])
      })
      jack.e.Y.est <- sapply(1:(length(U) - 1), function(k) {
        jack_Y[k] / get.NWK(jack_Y, jack_U, jack_U[k])
      })
      
      cov.jack <- mean(jack.e.X.est * jack.e.Y.est) - mean(jack.e.X.est) * mean(jack.e.Y.est)
      sig.jack.X <- mean(jack.e.X.est^2) - mean(jack.e.X.est)^2
      sig.jack.Y <- mean(jack.e.Y.est^2) - mean(jack.e.Y.est)^2
      
      rho.j[j] <- length(U) * rho.e.est - (length(U) - 1) * (cov.jack / sqrt(sig.jack.X * sig.jack.Y))
    }
    
    #########JACKKNIFE EMPIRICAL LIKELIHOOD#################
    
    ci.jel <- findci(rho.j)
    upper.jel[jj] <- ci.jel$Up
    lower.jel[jj] <- ci.jel$Low
    length.jel[jj] <- ci.jel$Up - ci.jel$Low
    coverage.jel[jj] <- (ci.jel$Low <= rho[ii]) & (rho[ii] <= ci.jel$Up)
  }
  
  # Store results for this rho
  results[ii, ] <- c(mean(lower.jel), mean(upper.jel), mean(length.jel), mean(coverage.jel))
}

# Display Results
data.frame(
  Rho = rho,
  Lower_Bound = results[, 1],
  Upper_Bound = results[, 2],
  CI_Length = results[, 3],
  Coverage_Probability = results[, 4]
)
