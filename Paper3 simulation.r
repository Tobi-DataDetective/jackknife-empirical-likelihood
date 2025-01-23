library(MASS)
library(emplik)
library(kedd)
library(progress)  # For progress bar

######### PARAMETERS ########################
n_values <- c(25, 50, 75, 100)  # Sample sizes
rho <- c(-0.9, -0.5, 0, 0.5, 0.9)  # Correlation coefficients
iterations <- 2000  # Number of iterations

########Define Matrix of Results ###################
results <- list()
colnames_results <- c("Lower Bound", "Upper Bound", "CI Length", "Coverage Probability")

########## FUNCTIONS #########################
get.NWK <- function(x, u, small.u) {
  bw <- density(u, kernel = "epanechnikov")$bw
  d <- outer(u, small.u, "-") / bw
  K <- (bw^(-1) * pmax(0, 0.75 * (1 - d^2)))
  result <- if (is.matrix(K) && all(colSums(K) > 0)) colSums(x * K) / colSums(K) else NA
  return(result)
}

findci <- function(x.vector) {
  cifunc <- function(ci.val, x.vector) {
    el.test(x.vector - ci.val, 0)
  }
  vectorMean <- mean(x.vector)
  L <- findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
  U <- findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
  return(list(Lower = L, Upper = U))
}

######### SIMULATE U FROM MULTIPLE DISTRIBUTIONS #########
generate_U <- function(n, distribution) {
  if (distribution == "uniform") {
    return(runif(n, 0, 1))
  } else if (distribution == "normal") {
    return(pnorm(rnorm(n, 0, 1)))
  } else if (distribution == "beta") {
    return(rbeta(n, 2, 5))
  } else if (distribution == "weibull") {
    return(rweibull(n, shape = 2, scale = 1))
  } else {
    stop("Unsupported distribution")
  }
}

######### SIMULATION #########################
# Initialize progress bar
total_steps <- length(n_values) * length(rho) * length(c("uniform", "normal", "beta", "weibull")) * iterations
pb <- progress_bar$new(
  format = "  Progress [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = total_steps, clear = FALSE, width = 60
)

for (n in n_values) {
  for (r in rho) {
    for (distribution in c("uniform", "normal", "beta", "weibull")) {
      
      lower.jel <- vector()
      upper.jel <- vector()
      coverage.jel <- vector()
      length.jel <- vector()
      
      for (iter in 1:iterations) {
        pb$tick()  # Update the progress bar
        
        ####### DATA GENERATION #####################
        sigma <- matrix(c(1, r, r, 1), nrow = 2)
        data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
        colnames(data) <- c("X", "Y")
        
        U <- generate_U(n, distribution)
        
        ######## DISTORTION FUNCTIONS ################
        if (distribution == "uniform") {
          ps_U <- 13 - U^2
          ph_U <- U^3 - 14
        } else if (distribution == "normal") {
          ps_U <- 0.5 * (U - mean(U))
          ph_U <- -U^2 + 0.25
        } else if (distribution == "beta") {
          ps_U <- U^3 - 0.5 * U
          ph_U <- -0.75 * U^2 + 0.1
        } else if (distribution == "weibull") {
          ps_U <- log(U + 0.5)
          ph_U <- U^1.5 - 0.25
        }
        
        ###### OBSERVED DATA #########################
        X_tilde <- data[, "X"] + ps_U
        Y_tilde <- data[, "Y"] + ph_U
        
        ####### CALCULATING ESTIMATORS #################
        e.X.est <- sapply(1:n, function(j) {
          result <- get.NWK(X_tilde, U, U[j])
          if (!is.na(result)) X_tilde[j] / result else NA
        })
        e.Y.est <- sapply(1:n, function(j) {
          result <- get.NWK(Y_tilde, U, U[j])
          if (!is.na(result)) Y_tilde[j] / result else NA
        })
        
        cov.e.est <- mean(e.X.est * e.Y.est, na.rm = TRUE) - mean(e.X.est, na.rm = TRUE) * mean(e.Y.est, na.rm = TRUE)
        sig.e.X <- mean(e.X.est^2, na.rm = TRUE) - mean(e.X.est, na.rm = TRUE)^2
        sig.e.Y <- mean(e.Y.est^2, na.rm = TRUE) - mean(e.Y.est, na.rm = TRUE)^2
        rho.e.est <- cov.e.est / sqrt(sig.e.X * sig.e.Y)
        
        ####### JACKKNIFING ###########################
        rho.j <- numeric(n)
        for (j in 1:n) {
          jack_X <- X_tilde[-j]
          jack_Y <- Y_tilde[-j]
          jack_U <- U[-j]
          
          jack.e.X.est <- sapply(1:(n - 1), function(k) {
            result <- get.NWK(jack_X, jack_U, jack_U[k])
            if (!is.na(result)) jack_X[k] / result else NA
          })
          jack.e.Y.est <- sapply(1:(n - 1), function(k) {
            result <- get.NWK(jack_Y, jack_U, jack_U[k])
            if (!is.na(result)) jack_Y[k] / result else NA
          })
          
          cov.jack <- mean(jack.e.X.est * jack.e.Y.est, na.rm = TRUE) - mean(jack.e.X.est, na.rm = TRUE) * mean(jack.e.Y.est, na.rm = TRUE)
          sig.jack.X <- mean(jack.e.X.est^2, na.rm = TRUE) - mean(jack.e.X.est, na.rm = TRUE)^2
          sig.jack.Y <- mean(jack.e.Y.est^2, na.rm = TRUE) - mean(jack.e.Y.est, na.rm = TRUE)^2
          rho.j[j] <- n * rho.e.est - (n - 1) * (cov.jack / sqrt(sig.jack.X * sig.jack.Y))
        }
        
        ######### EMPIRICAL LIKELIHOOD #################
        if (all(!is.na(rho.j) & is.finite(rho.j))) {
          ci.jel <- findci(rho.j)
        } else {
          ci.jel <- list(Upper = NA, Lower = NA)
        }
        
        upper.jel[iter] <- ci.jel$Upper
        lower.jel[iter] <- ci.jel$Lower
        length.jel[iter] <- if (!is.na(ci.jel$Upper) & !is.na(ci.jel$Lower)) ci.jel$Upper - ci.jel$Lower else NA
        coverage.jel[iter] <- if (!is.na(ci.jel$Lower) & !is.na(ci.jel$Upper)) (ci.jel$Lower <= r) & (r <= ci.jel$Upper) else NA
      }
      
      # Store results for this combination
      results[[paste0("n=", n, "_rho=", r, "_dist=", distribution)]] <- c(
        mean(lower.jel, na.rm = TRUE),
        mean(upper.jel, na.rm = TRUE),
        mean(length.jel, na.rm = TRUE),
        mean(coverage.jel, na.rm = TRUE)
      )
    }
  }
}

# Convert results to a tabular format
results_df <- do.call(rbind, lapply(names(results), function(key) {
  result <- results[[key]]
  data.frame(
    Scenario = key,
    Lower_Bound = round(result[1], 3),
    Upper_Bound = round(result[2], 3),
    CI_Length = round(result[3], 3),
    Coverage_Probability = round(result[4], 3)
  )
}))

# Display Results
print(results_df)



# 
# > library(MASS)
# > library(emplik)
# > library(kedd)
# > library(progress)  # For progress bar
# > 
#   > ######### PARAMETERS ########################
# > n_values <- c(25, 50, 75, 100)  # Sample sizes
# > rho <- c(-0.9, -0.5, 0, 0.5, 0.9)  # Correlation coefficients
# > iterations <- 2000  # Number of iterations
# > 
#   > ########Define Matrix of Results ###################
# > results <- list()
# > colnames_results <- c("Lower Bound", "Upper Bound", "CI Length", "Coverage Probability")
# > 
#   > ########## FUNCTIONS #########################
# > get.NWK <- function(x, u, small.u) {
#   +  bw <- density(u, kernel = "epanechnikov")$bw
#   +  d <- outer(u, small.u, "-") / bw
#   +  K <- (bw^(-1) * pmax(0, 0.75 * (1 - d^2)))
#   +  result <- if (is.matrix(K) && all(colSums(K) > 0)) colSums(x * K) / colSums(K) else NA
#   +  return(result)
#   + }
# > 
#   > findci <- function(x.vector) {
#     +  cifunc <- function(ci.val, x.vector) {
#       +    el.test(x.vector - ci.val, 0)
#       +  }
#     +  vectorMean <- mean(x.vector)
#     +  L <- findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Low
#     +  U <- findUL2(fun = cifunc, MLE = vectorMean, x.vector = x.vector)$Up
#     +  return(list(Lower = L, Upper = U))
#     + }
# > 
#   > ######### SIMULATE U FROM MULTIPLE DISTRIBUTIONS #########
# > generate_U <- function(n, distribution) {
#   +  if (distribution == "uniform") {
#     +    return(runif(n, 0, 1))
#     +  } else if (distribution == "normal") {
#       +    return(pnorm(rnorm(n, 0, 1)))
#       +  } else if (distribution == "beta") {
#         +    return(rbeta(n, 2, 5))
#         +  } else if (distribution == "weibull") {
#           +    return(rweibull(n, shape = 2, scale = 1))
#           +  } else {
#             +    stop("Unsupported distribution")
#             +  }
#   + }
# > 
#   > ######### SIMULATION #########################
# > # Initialize progress bar
#   > total_steps <- length(n_values) * length(rho) * length(c("uniform", "normal", "beta", "weibull")) * iterations
# > pb <- progress_bar$new(
#   +  format = "  Progress [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
#   +  total = total_steps, clear = FALSE, width = 60
#   + )
# > 
#   > for (n in n_values) {
#     +  for (r in rho) {
#       +    for (distribution in c("uniform", "normal", "beta", "weibull")) {
#         +      
#           +      lower.jel <- vector()
#           +      upper.jel <- vector()
#           +      coverage.jel <- vector()
#           +      length.jel <- vector()
#           +      
#             +      for (iter in 1:iterations) {
#               +        pb$tick()  # Update the progress bar
#               +        
#                 +        ####### DATA GENERATION #####################
#               +        sigma <- matrix(c(1, r, r, 1), nrow = 2)
#               +        data <- mvrnorm(n, mu = c(0, 0), Sigma = sigma)
#               +        colnames(data) <- c("X", "Y")
#               +        
#                 +        U <- generate_U(n, distribution)
#                 +        
#                   +        ######## DISTORTION FUNCTIONS ################
#                 +        if (distribution == "uniform") {
#                   +          ps_U <- 13 - U^2
#                   +          ph_U <- U^3 - 14
#                   +        } else if (distribution == "normal") {
#                     +          ps_U <- 0.5 * (U - mean(U))
#                     +          ph_U <- -U^2 + 0.25
#                     +        } else if (distribution == "beta") {
#                       +          ps_U <- U^3 - 0.5 * U
#                       +          ph_U <- -0.75 * U^2 + 0.1
#                       +        } else if (distribution == "weibull") {
#                         +          ps_U <- log(U + 0.5)
#                         +          ph_U <- U^1.5 - 0.25
#                         +        }
#                 +        
#                   +        ###### OBSERVED DATA #########################
#                 +        X_tilde <- data[, "X"] + ps_U
#                 +        Y_tilde <- data[, "Y"] + ph_U
#                 +        
#                   +        ####### CALCULATING ESTIMATORS #################
#                 +        e.X.est <- sapply(1:n, function(j) {
#                   +          result <- get.NWK(X_tilde, U, U[j])
#                   +          if (!is.na(result)) X_tilde[j] / result else NA
#                   +        })
#                 +        e.Y.est <- sapply(1:n, function(j) {
#                   +          result <- get.NWK(Y_tilde, U, U[j])
#                   +          if (!is.na(result)) Y_tilde[j] / result else NA
#                   +        })
#                 +        
#                   +        cov.e.est <- mean(e.X.est * e.Y.est, na.rm = TRUE) - mean(e.X.est, na.rm = TRUE) * mean(e.Y.est, na.rm = TRUE)
#                   +        sig.e.X <- mean(e.X.est^2, na.rm = TRUE) - mean(e.X.est, na.rm = TRUE)^2
#                   +        sig.e.Y <- mean(e.Y.est^2, na.rm = TRUE) - mean(e.Y.est, na.rm = TRUE)^2
#                   +        rho.e.est <- cov.e.est / sqrt(sig.e.X * sig.e.Y)
#                   +        
#                     +        ####### JACKKNIFING ###########################
#                   +        rho.j <- numeric(n)
#                   +        for (j in 1:n) {
#                     +          jack_X <- X_tilde[-j]
#                     +          jack_Y <- Y_tilde[-j]
#                     +          jack_U <- U[-j]
#                     +          
#                       +          jack.e.X.est <- sapply(1:(n - 1), function(k) {
#                         +            result <- get.NWK(jack_X, jack_U, jack_U[k])
#                         +            if (!is.na(result)) jack_X[k] / result else NA
#                         +          })
#                       +          jack.e.Y.est <- sapply(1:(n - 1), function(k) {
#                         +            result <- get.NWK(jack_Y, jack_U, jack_U[k])
#                         +            if (!is.na(result)) jack_Y[k] / result else NA
#                         +          })
#                       +          
#                         +          cov.jack <- mean(jack.e.X.est * jack.e.Y.est, na.rm = TRUE) - mean(jack.e.X.est, na.rm = TRUE) * mean(jack.e.Y.est, na.rm = TRUE)
#                         +          sig.jack.X <- mean(jack.e.X.est^2, na.rm = TRUE) - mean(jack.e.X.est, na.rm = TRUE)^2
#                         +          sig.jack.Y <- mean(jack.e.Y.est^2, na.rm = TRUE) - mean(jack.e.Y.est, na.rm = TRUE)^2
#                         +          rho.j[j] <- n * rho.e.est - (n - 1) * (cov.jack / sqrt(sig.jack.X * sig.jack.Y))
#                         +        }
#                   +        
#                     +        ######### EMPIRICAL LIKELIHOOD #################
#                   +        if (all(!is.na(rho.j) & is.finite(rho.j))) {
#                     +          ci.jel <- findci(rho.j)
#                     +        } else {
#                       +          ci.jel <- list(Upper = NA, Lower = NA)
#                       +        }
#                   +        
#                     +        upper.jel[iter] <- ci.jel$Upper
#                     +        lower.jel[iter] <- ci.jel$Lower
#                     +        length.jel[iter] <- if (!is.na(ci.jel$Upper) & !is.na(ci.jel$Lower)) ci.jel$Upper - ci.jel$Lower else NA
#                     +        coverage.jel[iter] <- if (!is.na(ci.jel$Lower) & !is.na(ci.jel$Upper)) (ci.jel$Lower <= r) & (r <= ci.jel$Upper) else NA
#                     +      }
#           +      
#             +      # Store results for this combination
#             +      results[[paste0("n=", n, "_rho=", r, "_dist=", distribution)]] <- c(
#               +        mean(lower.jel, na.rm = TRUE),
#               +        mean(upper.jel, na.rm = TRUE),
#               +        mean(length.jel, na.rm = TRUE),
#               +        mean(coverage.jel, na.rm = TRUE)
#               +      )
#             +    }
#       +  }
#     + }
# Progress [================] 100% | Elapsed:  5d | ETA:  0s
# > 
#   > # Convert results to a tabular format
#   > results_df <- do.call(rbind, lapply(names(results), function(key) {
#     +  result <- results[[key]]
#     +  data.frame(
#       +    Scenario = key,
#       +    Lower_Bound = round(result[1], 3),
#       +    Upper_Bound = round(result[2], 3),
#       +    CI_Length = round(result[3], 3),
#       +    Coverage_Probability = round(result[4], 3)
#       +  )
#     + }))
# > 
#   > # Display Results
#   > print(results_df)
# Scenario Lower_Bound Upper_Bound CI_Length Coverage_Probability
# 1   n=25_rho=-0.9_dist=uniform         NaN         NaN       NaN                  NaN
# 2    n=25_rho=-0.9_dist=normal         NaN         NaN       NaN                  NaN
# 3      n=25_rho=-0.9_dist=beta         NaN         NaN       NaN                  NaN
# 4   n=25_rho=-0.9_dist=weibull         NaN         NaN       NaN                  NaN
# 5   n=25_rho=-0.5_dist=uniform         NaN         NaN       NaN                  NaN
# 6    n=25_rho=-0.5_dist=normal         NaN         NaN       NaN                  NaN
# 7      n=25_rho=-0.5_dist=beta         NaN         NaN       NaN                  NaN
# 8   n=25_rho=-0.5_dist=weibull         NaN         NaN       NaN                  NaN
# 9      n=25_rho=0_dist=uniform         NaN         NaN       NaN                  NaN
# 10      n=25_rho=0_dist=normal         NaN         NaN       NaN                  NaN
# 11        n=25_rho=0_dist=beta         NaN         NaN       NaN                  NaN
# 12     n=25_rho=0_dist=weibull         NaN         NaN       NaN                  NaN
# 13   n=25_rho=0.5_dist=uniform         NaN         NaN       NaN                  NaN
# 14    n=25_rho=0.5_dist=normal         NaN         NaN       NaN                  NaN
# 15      n=25_rho=0.5_dist=beta         NaN         NaN       NaN                  NaN
# 16   n=25_rho=0.5_dist=weibull         NaN         NaN       NaN                  NaN
# 17   n=25_rho=0.9_dist=uniform         NaN         NaN       NaN                  NaN
# 18    n=25_rho=0.9_dist=normal         NaN         NaN       NaN                  NaN
# 19      n=25_rho=0.9_dist=beta         NaN         NaN       NaN                  NaN
# 20   n=25_rho=0.9_dist=weibull         NaN         NaN       NaN                  NaN
# 21  n=50_rho=-0.9_dist=uniform         NaN         NaN       NaN                  NaN
# 22   n=50_rho=-0.9_dist=normal         NaN         NaN       NaN                  NaN
# 23     n=50_rho=-0.9_dist=beta         NaN         NaN       NaN                  NaN
# 24  n=50_rho=-0.9_dist=weibull         NaN         NaN       NaN                  NaN
# 25  n=50_rho=-0.5_dist=uniform         NaN         NaN       NaN                  NaN
# 26   n=50_rho=-0.5_dist=normal         NaN         NaN       NaN                  NaN
# 27     n=50_rho=-0.5_dist=beta         NaN         NaN       NaN                  NaN
# 28  n=50_rho=-0.5_dist=weibull         NaN         NaN       NaN                  NaN
# 29     n=50_rho=0_dist=uniform         NaN         NaN       NaN                  NaN
# 30      n=50_rho=0_dist=normal         NaN         NaN       NaN                  NaN
# 31        n=50_rho=0_dist=beta         NaN         NaN       NaN                  NaN
# 32     n=50_rho=0_dist=weibull         NaN         NaN       NaN                  NaN
# 33   n=50_rho=0.5_dist=uniform         NaN         NaN       NaN                  NaN
# 34    n=50_rho=0.5_dist=normal         NaN         NaN       NaN                  NaN
# 35      n=50_rho=0.5_dist=beta         NaN         NaN       NaN                  NaN
# 36   n=50_rho=0.5_dist=weibull         NaN         NaN       NaN                  NaN
# 37   n=50_rho=0.9_dist=uniform         NaN         NaN       NaN                  NaN
# 38    n=50_rho=0.9_dist=normal         NaN         NaN       NaN                  NaN
# 39      n=50_rho=0.9_dist=beta         NaN         NaN       NaN                  NaN
# 40   n=50_rho=0.9_dist=weibull         NaN         NaN       NaN                  NaN
# 41  n=75_rho=-0.9_dist=uniform         NaN         NaN       NaN                  NaN
# 42   n=75_rho=-0.9_dist=normal         NaN         NaN       NaN                  NaN
# 43     n=75_rho=-0.9_dist=beta         NaN         NaN       NaN                  NaN
# 44  n=75_rho=-0.9_dist=weibull         NaN         NaN       NaN                  NaN
# 45  n=75_rho=-0.5_dist=uniform         NaN         NaN       NaN                  NaN
# 46   n=75_rho=-0.5_dist=normal         NaN         NaN       NaN                  NaN
# 47     n=75_rho=-0.5_dist=beta         NaN         NaN       NaN                  NaN
# 48  n=75_rho=-0.5_dist=weibull         NaN         NaN       NaN                  NaN
# 49     n=75_rho=0_dist=uniform         NaN         NaN       NaN                  NaN
# 50      n=75_rho=0_dist=normal         NaN         NaN       NaN                  NaN
# 51        n=75_rho=0_dist=beta         NaN         NaN       NaN                  NaN
# 52     n=75_rho=0_dist=weibull         NaN         NaN       NaN                  NaN
# 53   n=75_rho=0.5_dist=uniform         NaN         NaN       NaN                  NaN
# 54    n=75_rho=0.5_dist=normal         NaN         NaN       NaN                  NaN
# 55      n=75_rho=0.5_dist=beta         NaN         NaN       NaN                  NaN
# 56   n=75_rho=0.5_dist=weibull         NaN         NaN       NaN                  NaN
# 57   n=75_rho=0.9_dist=uniform         NaN         NaN       NaN                  NaN
# 58    n=75_rho=0.9_dist=normal         NaN         NaN       NaN                  NaN
# 59      n=75_rho=0.9_dist=beta         NaN         NaN       NaN                  NaN
# 60   n=75_rho=0.9_dist=weibull         NaN         NaN       NaN                  NaN
# 61 n=100_rho=-0.9_dist=uniform         NaN         NaN       NaN                  NaN
# 62  n=100_rho=-0.9_dist=normal         NaN         NaN       NaN                  NaN
# 63    n=100_rho=-0.9_dist=beta         NaN         NaN       NaN                  NaN
# 64 n=100_rho=-0.9_dist=weibull         NaN         NaN       NaN                  NaN
# 65 n=100_rho=-0.5_dist=uniform         NaN         NaN       NaN                  NaN
# 66  n=100_rho=-0.5_dist=normal         NaN         NaN       NaN                  NaN
# 67    n=100_rho=-0.5_dist=beta         NaN         NaN       NaN                  NaN
# 68 n=100_rho=-0.5_dist=weibull         NaN         NaN       NaN                  NaN
# 69    n=100_rho=0_dist=uniform         NaN         NaN       NaN                  NaN
# 70     n=100_rho=0_dist=normal         NaN         NaN       NaN                  NaN
# 71       n=100_rho=0_dist=beta         NaN         NaN       NaN                  NaN
# 72    n=100_rho=0_dist=weibull         NaN         NaN       NaN                  NaN
# 73  n=100_rho=0.5_dist=uniform         NaN         NaN       NaN                  NaN
# 74   n=100_rho=0.5_dist=normal         NaN         NaN       NaN                  NaN
# 75     n=100_rho=0.5_dist=beta         NaN         NaN       NaN                  NaN
# 76  n=100_rho=0.5_dist=weibull         NaN         NaN       NaN                  NaN
# 77  n=100_rho=0.9_dist=uniform         NaN         NaN       NaN                  NaN
# 78   n=100_rho=0.9_dist=normal         NaN         NaN       NaN                  NaN
# 79     n=100_rho=0.9_dist=beta         NaN         NaN       NaN                  NaN
# 80  n=100_rho=0.9_dist=weibull         NaN         NaN       NaN                  NaN
# > results
# $`n=25_rho=-0.9_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=-0.9_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=-0.9_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=-0.9_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=-0.5_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=-0.5_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=-0.5_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=-0.5_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=0_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=0_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=0_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=0_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=0.5_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=0.5_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=0.5_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=0.5_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=0.9_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=0.9_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=0.9_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=25_rho=0.9_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=-0.9_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=-0.9_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=-0.9_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=-0.9_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=-0.5_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=-0.5_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=-0.5_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=-0.5_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=0_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=0_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=0_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=0_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=0.5_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=0.5_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=0.5_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=0.5_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=0.9_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=0.9_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=0.9_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=50_rho=0.9_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=-0.9_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=-0.9_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=-0.9_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=-0.9_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=-0.5_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=-0.5_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=-0.5_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=-0.5_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=0_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=0_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=0_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=0_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=0.5_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=0.5_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=0.5_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=0.5_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=0.9_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=0.9_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=0.9_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=75_rho=0.9_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=-0.9_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=-0.9_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=-0.9_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=-0.9_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=-0.5_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=-0.5_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=-0.5_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=-0.5_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=0_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=0_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=0_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=0_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=0.5_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=0.5_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=0.5_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=0.5_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=0.9_dist=uniform`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=0.9_dist=normal`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=0.9_dist=beta`
# [1] NaN NaN NaN NaN
# 
# $`n=100_rho=0.9_dist=weibull`
# [1] NaN NaN NaN NaN
# 
# 

