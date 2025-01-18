library(MASS)
library(emplik)
library(kedd)

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
for (n in n_values) {
  for (r in rho) {
    for (distribution in c("uniform", "normal", "beta", "weibull")) {
      
      lower.jel <- vector()
      upper.jel <- vector()
      coverage.jel <- vector()
      length.jel <- vector()
      
      for (iter in 1:iterations) {
        
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
          print(list(rho.j = rho.j, r = r, n = n, iter = iter))  # Debug output
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

# Optionally save the results to a CSV file
# write.csv(results_df, "paper3_simulation_results.csv", row.names = FALSE)
