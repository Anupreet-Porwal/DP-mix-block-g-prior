# Load necessary libraries
library(invgamma)       # For inverse gamma distribution functions
library(WoodburyMatrix) # For efficient matrix operations on Woodbury matrices
library(mvtnorm)        # For multivariate normal distributions
library(truncnorm)      # For truncated normal distribution functions
library(extraDistr)     # For additional distribution functions
library(dirichletprocess) # For Dirichlet process mixture modeling

# Function to resample elements of a vector
# Parameters:
#   x: the vector to be resampled
#   ...: additional arguments passed to sample.int
# Returns:
#   A resampled version of vector x
resample <- function(x, ...) x[sample.int(length(x), ...)]

# Function to transform a vector 'a' into a lexicographic order starting from 1
# Parameters:
#   a: numeric vector to be reordered
# Returns:
#   Vector 'a' with unique elements re-indexed from 1 in lexicographic order
get_lexi_order2 <- function(a) {
  # Return 'a' if empty
  if (length(a) == 0) { return(a) }
  
  # Get unique values of 'a'
  y <- unique(a)
  
  # Check if 'y' is already lexicographically ordered from 1
  if (any(diff(y) > 1) | min(y) != 1 | y[1] != 1) {
    # Re-index elements of 'a' to start from 1
    af1 <- factor(a, levels = y)       # Create factor with unique values of 'a'
    levels(af1) <- seq(1, length(y))   # Set levels to 1 through the number of unique values
    return(as.numeric(af1))            # Return as numeric vector
  } else {
    return(a)                          # Return 'a' if already ordered
  }
}


# Function to calculate the log-likelihood of cluster proportions
# Parameters:
#   k: cluster index
#   b: coefficient vector
#   xtx: cross-product matrix of predictors
#   g: vector of gamma parameters for shrinkage
#   gam: binary vector indicating variable inclusion (1) or exclusion (0)
#   sigma2: variance parameter
#   g_K: vector of cluster-specific shrinkage parameters
#   r: current variable index for likelihood calculation
#   prior.corr: logical, whether to consider prior correlation
#   tau2: additional variance scaling parameter
# Returns:
#   Calculated log-likelihood of the cluster proportion
cluster_prop_llikelihood <- function(k, b, xtx, g, gam, sigma2, g_K, r, prior.corr, tau2) {
  var_in_mod <- which(gam == 1)          # Indices of variables included in the model
  ind.var <- var_in_mod[r]               # Select the r-th included variable
  g[ind.var] <- g_K[k]                   # Assign cluster-specific shrinkage parameter
  ggam <- g[as.logical(gam)]             # Shrinkage parameters for included variables
  
  # Calculate transformed coefficient vector for likelihood calculation
  b_sqrtg <- t(b) %*% diag(1 / sqrt(ggam), nrow = length(ggam), ncol = length(ggam))
  
  # Calculate individual likelihood component based on prior correlation setting
  if (prior.corr == TRUE) {
    d_ind <- 2 * (xtx[r, -r] %*% b_sqrtg[1, -r]) * b[r] / sqrt(ggam[r]) +
      xtx[r, r] * b[r]^2 / ggam[r]
  } else {
    d_ind <- b[r]^2 / ggam[r]
  }
  
  # Calculate final log-likelihood
  prop_llik <- -0.5 * log(ggam[r]) - d_ind / (2 * sigma2 * tau2)
  
  return(prop_llik)
}


# Simplified version of the cluster proportion log-likelihood calculation
# Parameters:
#   k: cluster index
#   b: coefficient vector
#   xtx: cross-product matrix of predictors
#   sigma2: variance parameter
#   g_K: vector of cluster-specific shrinkage parameters
#   grp_idx.gam: grouping index for gamma parameters
#   r: current variable index for likelihood calculation
#   prior.corr: logical, whether to consider prior correlation
#   tau2: additional variance scaling parameter
# Returns:
#   Calculated simplified log-likelihood of the cluster proportion
cluster_prop_llikelihood_simplified <- function(k, b, xtx, sigma2, g_K, grp_idx.gam, r, prior.corr, tau2) {
  grp_idx.gam[r] <- k                    # Assign cluster index to current variable
  ggam <- g_K[grp_idx.gam]               # Shrinkage parameters for group index
  
  # Calculate transformed coefficient vector for likelihood calculation
  b_sqrtg <- t(b) %*% diag(1 / sqrt(ggam), nrow = length(ggam), ncol = length(ggam))
  sqrtg_diag <- diag(sqrt(ggam), nrow = length(ggam), ncol = length(ggam))
  
  # Calculate individual likelihood component based on prior correlation setting
  if (prior.corr == TRUE) {
    d_ind <- 2 * (xtx[r, -r] %*% b_sqrtg[1, -r]) * b[r] / sqrt(ggam[r]) +
      xtx[r, r] * b[r]^2 / ggam[r]
  } else {
    d_ind <- b[r]^2 / ggam[r]
  }
  
  # Calculate final log-likelihood
  prop_llik <- -0.5 * log(ggam[r]) - d_ind / (2 * sigma2 * tau2)
  
  return(prop_llik)
}

# Function to compute posterior probability of cluster proportion
# Parameters:
#   k: cluster index
#   b: coefficient vector
#   xtx: cross-product matrix of predictors
#   sigma2: variance parameter
#   g_K: vector of cluster-specific shrinkage parameters
#   grp_idx.gam: grouping index for gamma parameters
#   r: current variable index for likelihood calculation
#   prior.corr: logical, whether to consider prior correlation
#   tau2: additional variance scaling parameter
# Returns:
#   Log-posterior probability of cluster proportion
cluster_prop_post <- function(k, b, xtx, sigma2, g_K, grp_idx.gam, r, prior.corr, tau2) {
  # Calculate posterior as sum of cluster count log and simplified log-likelihood
  return(log(sum(grp_idx.gam == k)) + cluster_prop_llikelihood_simplified(k, b, xtx, sigma2, g_K, grp_idx.gam, r, prior.corr, tau2))
}


# Taken from: https://github.com/IvanUkhov/blog/blob/main/_scripts/2021-01-25-dirichlet-process/common.R
# blog post: https://blog.ivanukhov.com/2021/01/25/dirichlet-process.html

# Function to count the number of subjects in each group
# Parameters:
#   m: total number of groups
#   k: vector of group labels for each subject
# Returns:
#   Vector 'n' with the count of subjects in each group
count_subjects <- function(m, k) {
  n <- rep(0, m)                            # Initialize count vector with zeroes
  if (length(k) == 0) {                     # If 'k' is empty, return 'n' as is
    return(n)
  }
  k <- as.data.frame(table(k))              # Create a frequency table of 'k'
  n[as.numeric(levels(k[, 1]))[k[, 1]]] <- k[, 2] # Assign counts to corresponding indices
  return(n)
}


# Function to sample stick-breaking weights for a Dirichlet process
# Parameters:
#   l: number of breaks
#   alpha, beta: parameters of the Beta distribution
#   log: logical, if TRUE returns log of stick-breaking weights
# Returns:
#   List containing 'p' (stick-breaking weights) and 'q' (Beta samples)
stick_break <- function(l, alpha, beta, log = TRUE) {
  q <- rbeta(l, shape1 = alpha, shape2 = beta)     # Sample from Beta distribution
  q <- c(head(q, -1), 1)                           # Ensure last component is 1
  
  # Compute stick-breaking weights in log or standard scale
  if (log == TRUE) {
    p <- log(q) + c(0, cumsum(log(1 - head(q, -1)))) # Log-scale weights
  } else {
    p <- q * c(1, cumprod(1 - head(q, -1)))          # Standard scale weights
  }
  
  return(list(p = p, q = q))
}


# Function to sample values from the prior distribution of lambda
# Parameters:
#   l: number of lambda values to sample
#   alpha0, beta0: shape and rate parameters for the Gamma distribution
# Returns:
#   Vector of lambda samples from the Gamma prior
sample_Plambda_prior <- function(l, alpha0 = 1, beta0 = 1) {
  rgamma(l, alpha0, beta0)  # Sample from Gamma(alpha0, beta0) distribution
}

# Function to sample values from the posterior distribution of lambda
# Parameters:
#   l: number of lambda values to sample
#   q: vector of stick-breaking proportions
#   alpha0, beta0: shape and rate parameters for the Gamma distribution
# Returns:
#   Vector of lambda samples from the Gamma posterior
sample_Plambda_posterior <- function(l, q, alpha0 = 1, beta0 = 1) {
  sample_Plambda_prior(
    l = l,
    alpha0 = alpha0 + length(q) - 1,               # Update alpha with prior and observed q length
    beta0 = beta0 - sum(log(head(q, -1)))          # Update beta with the log-sum of q values
  )
}



# Taken from here :  https://github.com/willtownes/mit6882/blob/master/will/speed_tests.Rmd
# Some discussion here:  http://www.statsathome.com/2018/10/19/sampling-from-multivariate-normal-precision-and-covariance-parameterizations/

# Function to generate random samples based on custom transformation of a multivariate normal distribution
# Parameters:
#   n: number of samples
#   theta: mean vector
#   Lambda: covariance matrix
#   D: dimension of theta (defaults to length of theta)
# Returns:
#   Transformed matrix of samples
custom4 <- function(n, theta, Lambda, D = length(theta)) {
  Q <- chol(Lambda)                      # Cholesky decomposition of Lambda
  Z <- matrix(rnorm(n * D), nrow = D, ncol = n)  # Generate standard normal random matrix
  # Transformation: Q^{-1} * Z + Lambda^{-1} * theta
  backsolve(Q, Z + drop(backsolve(Q, theta, transpose = TRUE)))
}

# Sample from a hyper-g prior 
rhyperg <- function(n,a){
  h <- rbeta(n,1,a/2-1)
  
  while(any(h/(1-h)==Inf)){
    h <- rbeta(n,1,a/2-1)
  }
  return(h/(1-h))
}

# Sample from a hyper-g-n prior 
rhypergn <- function(x,a,n){
  h <- rbeta(x,1,a/2-1)
  return(n*h/(1-h))
}


# Function to sample from a truncated Gamma distribution
# Parameters:
#   n: number of samples
#   shape, rate: parameters for the Gamma distribution
#   truncation: truncation point
# Returns:
#   Vector of samples truncated at the specified point
rtgamma <- function(n, shape, rate, truncation) {
  # Step 1: Find probability corresponding to truncation point
  p_trunc <- pgamma(truncation, shape = shape, rate = rate)
  
  # Step 2: Draw uniform samples within the truncation probability range
  u_truncgam <- runif(n, 0, p_trunc)
  
  # Step 3: Transform using the quantile function for Gamma distribution
  t <- qgamma(u_truncgam, shape = shape, rate = rate)
  
  # Ensure no zero or infinite values in output
  while (t == 0 | t == Inf) {
    u_truncgam <- runif(n, 0, p_trunc)
    t <- qgamma(u_truncgam, shape = shape, rate = rate)
  }
  
  return(t)
}

# Function to propose a new gamma vector by modifying its entries
# Parameters:
#   gam: binary vector representing the current gamma state
#   d: maximum number of entries to change (default is 4)
# Returns:
#   Modified gamma vector
proposal.gamma <- function(gam, d = 4) {
  p <- length(gam)                            # Length of the gamma vector
  prob1 <- c(1, 0, 0, 0)                      # Probability for number of entries to modify; Currently set to add/remove only one variable at a time 
  prob2 <- c(0.7, 0.3)                        # Probability for choice of modification type
  
  # Step 1: Randomly choose between two modification strategies
  s <- resample(1:2, 1, prob = prob2)
  
  # Strategy 1: Flip a subset of entries in gamma
  if (s == 1) {
    # Adjust d to be within the length of gamma
    if (p < d) {
      d <- p
      prob1 <- prob1[1:d] / sum(prob1[1:d])   # Normalize prob1 if d is modified
    }
    
    # Randomly select the number of entries to flip and their indices
    d0 <- resample(1:d, 1, prob = prob1)
    ind <- resample(1:p, d0)
    gam[ind] =!gam[ind]                     # Flip selected entries
    
    # Strategy 2: Swap an entry from active set with one from non-active set
  } else if (s == 2) {
    if (all(gam == 1)) {                      # Case: all entries are active (1)
      ind <- resample(1:p, 1)
      gam[ind] <- 0                           # Randomly set one entry to inactive (0)
      
    } else if (all(gam == 0)) {               # Case: all entries are inactive (0)
      ind <- resample(1:p, 1)
      gam[ind] <- 1                           # Randomly set one entry to active (1)
      
    } else {                                  # Case: mixed active and inactive entries
      gam1 <- which(gam == 1)                 # Indices of active entries
      gam0 <- which(gam == 0)                 # Indices of inactive entries
      ind1 <- resample(gam1, 1)               # Choose one active entry
      ind0 <- resample(gam0, 1)               # Choose one inactive entry
      gam[ind1] <- 0                          # Swap: set chosen active entry to inactive
      gam[ind0] <- 1                          # Swap: set chosen inactive entry to active
    }
  }
  
  return(gam)
}



# Function to compute the log marginal likelihood for a Bayesian linear model
# Parameters:
#   x: predictor matrix
#   y: response vector
#   g_K: vector of group-level scale parameters
#   grp_idx: vector of group indices for predictors
#   gam: binary vector indicating active predictors
#   hyper.prior: prior type, e.g., "beta-prime-MG"
#   prior.corr: logical, if TRUE uses correlation-based prior
#   tau2: prior scale parameter
# Returns:
#   A list with the log marginal likelihood and intermediate matrices
log_marginal <- function(x, y, g_K, grp_idx, gam, hyper.prior, prior.corr, tau2) {
  
  n <- length(y)             # Number of observations
  p <- ncol(x)               # Number of predictors
  
  # Subset x and grp_idx based on active predictors (where gam == 1)
  x <- x[, as.logical(gam)]
  grp_idx.gam <- grp_idx[as.logical(gam)]
  g <- g_K[grp_idx.gam]      # Scale parameters for active groups
  pgam <- sum(gam)           # Number of active predictors
  
  if (length(g) == 0) {
    # Case: No active predictors, use uninformative prior calculation
    num1 <- sum((y - mean(y))^2)
    logmg <- lgamma((n - 1) / 2) - 0.5 * log(n) - 
      0.5 * (n - 1) * log(pi) - 0.5 * (n - 1) * log(num1)
  } else {
    # Case: Active predictors, calculate using Woodbury matrix inversion
    xtx <- t(x) %*% x              # Cross-product of active predictors
    xty <- t(x) %*% y              # Predictor-response cross-product
    
    if (hyper.prior == "beta-prime-MG") {
      # For beta-prime hyperprior, rescale x using group-level scales
      ghalf <- diag(sqrt(g), nrow = pgam, ncol = pgam)
      xghalf <- x %*% ghalf
      
      # Construct Woodbury matrix with correlation adjustment if prior.corr is TRUE
      if (prior.corr == TRUE) {
        mat <- WoodburyMatrix(A = diag(n), X = sqrt(tau2) * xghalf, B = xtx)
      } else {
        mat <- WoodburyMatrix(A = diag(n), X = sqrt(tau2) * xghalf, B = diag(pgam))
      }
    } else {
      # Alternative prior type: use inverse scaling of xtx by g and tau2
      ghalf.inv <- diag(1 / sqrt(g), nrow = pgam, ncol = pgam)
      
      # Define inverse of scaled cross-product matrix
      if (prior.corr == TRUE) {
        gxtxg <- ghalf.inv %*% xtx %*% ghalf.inv  
      } else {
        gxtxg <- diag(1 / g, nrow = pgam, ncol = pgam)
      }
      
      gxtxg <- gxtxg / tau2
      gxtxg <- (gxtxg + t(gxtxg)) / 2  # Symmetrize matrix for numerical stability
      mat <- WoodburyMatrix(A = diag(n), X = x, B = gxtxg)
    }
    
    # Compute residual sum of squares with modified Woodbury matrix
    num1 <- as.numeric(y %*% solve(mat) %*% y - n * mean(y)^2)
    den <- as.numeric(determinant(mat, logarithm = TRUE)$modulus)
    
    logmg <- lgamma((n - 1) / 2) - 0.5 * log(n) - 
      0.5 * (n - 1) * log(pi) - 0.5 * (n - 1) * log(num1) - 0.5 * den
  }
  
  # Return results as a list including the log marginal likelihood and intermediate matrices
  res <- list("logmg" = logmg, "xtx" = xtx, "xty" = xty, 
              "ggam" = g, "g_K" = g_K, "grp_idx" = grp_idx, 
              "sig2rate" = num1)
  return(res)
}


#### Algorithm 1 of Liu et al 2013 ####

gamma_rejection_sampler <- function(alpha, gam_param, truncation = NULL) {
  
  # Calculate constant C based on gam_param
  if (gam_param != 0) {
    logC <- log(abs(gam_param)) - 0.5 * log(alpha)
    C <- sign(gam_param) * exp(logC)
  } else {
    C <- 0
  }
  
  c <- 0  # Counter for iterations
  
  # If truncation is not provided
  if (is.null(truncation)) {
    
    # Case when gam_param is less than or equal to 0
    if (gam_param <= 0) {
      logdelta0 <- log(4 * alpha) - 2 * log(sqrt(gam_param^2 + 4 * alpha) + abs(gam_param))
      delta0 <- exp(logdelta0)  # Calculate delta0
      
      # Sample from the gamma distribution
      t <- rgamma(1, shape = alpha, rate = delta0)
      u <- runif(1)  # Generate uniform random number
      
      # Acceptance-rejection loop
      while (log(u) > alpha - t * (1 - delta0) - 2 * sqrt(t) * gam_param - alpha / delta0) {
        t <- rgamma(1, shape = alpha, rate = delta0)
        u <- runif(1)
      }
      
    } else {  # Case when gam_param is greater than 0
      
      if (alpha > 1) {
        # Calculate upsilon and related values for sampling
        upsilon <- (-gam_param + sqrt(gam_param^2 + 4 * (alpha - 1)))^2 / 4
        xst = ((par[1] - upsilon) / par[2])^2
        logMst = (par[1] - upsilon) * log(xst) - 2 * par[2] * sqrt(xst)
        
        # Sample from the gamma distribution
        t <- rgamma(1, upsilon, 1)
        u <- runif(1)
        
        # Acceptance-rejection loop
        while (log(u) > (alpha - upsilon) * log(t) - 2 * gam_param * sqrt(t) - logMst) {
          t <- rgamma(1, upsilon, 1)
          u <- runif(1)
        }
        
      } else {  # Case when alpha is less than or equal to 1
        t <- rgamma(1, alpha, 1)
        u <- runif(1) 
        
        # Acceptance-rejection loop
        while (log(u) > -2 * gam_param * sqrt(t)) {
          t <- rgamma(1, alpha, 1)
          u <- runif(1)
        }
      }
      
    }
    
  } else {  # If truncation is provided
    if (gam_param <= 0) {
      logdelta0 <- log(4 * alpha) - 2 * log(sqrt(gam_param^2 + 4 * alpha) + abs(gam_param))
      delta0 <- exp(logdelta0)
      
      # Sample from the truncated gamma distribution
      t <- rtgamma(1, shape = alpha, rate = delta0, truncation = truncation)
      u <- runif(1) 
      
      # Acceptance-rejection loop
      while (log(u) > alpha - t * (1 - delta0) - 2 * sqrt(t) * gam_param - alpha / delta0) {
        c <- c + 1
        if (c > 1005) {
          # Debugging information before stopping
          print(paste("logu:", log(u)))
          print(paste("compared expression:", alpha - t * (1 - delta0) - 2 * sqrt(t) * gam_param - alpha / delta0))
          print(paste("last gamma val:", t))
          print(paste("alpha:", alpha))
          print(paste("gam_param:", gam_param))
          print(paste("delta0:", delta0))
          print(paste("truncation:", truncation))
          print("------")
          
          stop("Stuck in algo 1, gam < 0")
        }
        t <- rtgamma(1, shape = alpha, rate = delta0, truncation = truncation)
        u <- runif(1)
      }
      
    } else {  # Case when gam_param is greater than 0
      if (alpha > 1) {
        par <- c(alpha, gam_param)
        upsilon <- (-gam_param + sqrt(gam_param^2 + 4 * (alpha - 1)))^2 / 4
        xst = min(((par[1] - upsilon) / par[2])^2, truncation)
        logMst = (par[1] - upsilon) * log(xst) - 2 * par[2] * sqrt(xst)
        
        # Sample from the truncated gamma distribution
        t <- rtgamma(1, shape = upsilon, rate = 1, truncation = truncation)
        u <- runif(1) 
        
        # Acceptance-rejection loop
        while (log(u) > (alpha - upsilon) * log(t) - 2 * gam_param * sqrt(t) - logMst) {
          c <- c + 1
          t <- rtgamma(1, shape = upsilon, rate = 1, truncation = truncation)
          u <- runif(1)
        }
        
      } else {  # Case when alpha is less than or equal to 1
        t <- rtgamma(1, shape = alpha, rate = 1, truncation = truncation)
        u <- runif(1) 
        
        # Acceptance-rejection loop
        while (log(u) > -2 * gam_param * sqrt(t)) {
          c <- c + 1
          t <- rtgamma(1, shape = alpha, rate = 1, truncation = truncation)
          u <- runif(1)
        }
      }
    }
  }
  
  return(t)  # Return the sampled value
}


#### Algorithm 2 of Liu et al 2013 ####
normal_rejection_sampler <- function(alpha, gam_param, truncation = NULL) {
  c <- 0  # Counter for iterations
  
  # Check if alpha is less than 1/2
  if (alpha < 1/2) {
    stop("For normal RS, alpha needs to be greater than 1/2")
  }
  
  # Calculate the maximum value of m based on alpha and gam_param
  if (alpha == 1/2) {
    m <- max(-gam_param, 0)  # When alpha is exactly 1/2
  } else {
    logm <- log(2 * alpha - 1) - log(gam_param + sqrt(gam_param^2 + 4 * alpha - 2))
    m <- exp(logm)  # Calculate m using the exponential of logm
  }
  
  # If no truncation is provided
  if (is.null(truncation)) {
    # Sample from truncated normal distribution
    x <- rtruncnorm(1, mean = m, sd = 1/sqrt(2), a = 0)
    u <- runif(1)  # Generate uniform random number
    
    # Acceptance-rejection loop
    while ((2 * alpha - 1) * log(m) + log(u) > (2 * alpha - 1) * log(x) - 2 * (m + gam_param) * (x - m)) {
      x <- rtruncnorm(1, mean = m, sd = 1/sqrt(2), a = 0)  # Resample x
      u <- runif(1)  # Resample u
    }
    
  } else {  # If truncation is provided
    # Sample from truncated normal distribution with upper bound
    x <- rtruncnorm(1, mean = m, sd = 1/sqrt(2), a = 0, b = sqrt(truncation))
    u <- runif(1)
    
    # Acceptance-rejection loop
    while ((2 * alpha - 1) * log(m) + log(u) > (2 * alpha - 1) * log(x) - 2 * (m + gam_param) * (x - m)) {
      c <- c + 1  # Increment the counter
      
      # Uncomment the following lines for debugging if needed
      if (c > 1000) {
       print(paste("lhs:", (2 * alpha - 1) * log(m) + log(u)))
       print(paste("rhs:", (2 * alpha - 1) * log(x) - 2 * (m + gam_param) * (x - m)))
       print(paste("last val:", x))
       print(paste("alpha:", alpha))
       print(paste("gam_param:", gam_param))
       print(paste("m:", m))
       print(paste("truncation:", truncation))
       print("------")
    
       if (c > 1005) {
         stop("Stuck in algo 2")
       }
     }
      
      # Resample x and u
      x <- rtruncnorm(1, mean = m, sd = 1/sqrt(2), a = 0, b = sqrt(truncation))
      u <- runif(1)
    }
  }
  
  return(x^2)  # Return the squared value of x
}


# Algorithm 3 of Liu et al 2013 ####
gamma_rs_sqrt_scale <- function(alpha, gam_param, truncation = NULL) {
  # Check if alpha is less than or equal to zero
  if (alpha < 0) {
    stop("For gamma RS on sqrt scale, alpha > 0")
  }
  
  # Calculate parameters for the gamma distribution
  delta1 <- gam_param + sqrt(gam_param^2 + 4 * alpha)  # Adjusted scale parameter
  r <- 2 * alpha  # Shape parameter for gamma distribution
  c <- 0  # Counter for iterations in rejection sampling
  
  # If no truncation is provided
  if (is.null(truncation)) {
    # Sample from the gamma distribution
    x <- rgamma(1, r, delta1)
    u <- runif(1)  # Generate a uniform random number
    
    # Acceptance-rejection loop
    while (log(u) > -(x - (delta1 / 2 - gam_param))^2) {
      x <- rgamma(1, r, delta1)  # Resample x
      u <- runif(1)  # Resample u
    }
  } else {  # If truncation is provided
    # Sample from the truncated gamma distribution
    x <- rtgamma(1, shape = r, rate = delta1, truncation = sqrt(truncation))
    u <- runif(1)
    
    # Acceptance-rejection loop
    while (log(u) > -(x - (delta1 / 2 - gam_param))^2) {
      c <- c + 1  # Increment the counter
      
      # Uncomment the following lines for debugging if needed
      if (c > 1000) {
        print(paste("logu:", log(u)))
        print(paste("rhs:", -(x - (delta1 / 2 - gam_param))^2))
        print(paste("last val:", x))
        print(paste("alpha:", alpha))
        print(paste("gam_param:", gam_param))
        print(paste("delta1:", delta1))
        print(paste("truncation:", truncation))
        print("------")
        
        # Stop the algorithm if it exceeds a certain number of iterations
        if (c > 1005) {
          stop("Stuck in algo 3")
        }
      }
      
      # Resample x and u
      x <- rtgamma(1, shape = r, rate = delta1, truncation = sqrt(truncation))
      u <- runif(1)
    }
  }
  
  return(x^2)  # Return the square of the sampled value
}

#### Optimal Rejection sampler for extended gamma distribution (Liu el al 2013) ####
ext_gamma_sampler <- function(alpha, gam_param, truncation = NULL) {
  # Check if gam_param is not zero
  if (gam_param != 0) {
    # Calculate the constant C for gamma sampling
    logC <- log(abs(gam_param)) - 0.5 * log(alpha)
    C <- sign(gam_param) * exp(logC)  # Determine the sign of gam_param
  } else if (gam_param == 0) {
    C <- 0  # Set C to 0 if gam_param is zero
  }
  
  # Check if alpha is NaN (not a number)
  if (is.nan(alpha)) {
    print("alpha is NaN")
  }
  
  # Choose the sampling method based on the value of C
  if (C <= -0.7) {
    # Use normal rejection sampling for C <= -0.7
    z <- normal_rejection_sampler(alpha, gam_param, truncation)
  } else if (-0.7 < C & C < 0.7) {
    # Use gamma rejection sampling for -0.7 < C < 0.7
    z <- gamma_rejection_sampler(alpha, gam_param, truncation)
  } else if (C >= 0.7) {
    # Use gamma sampling on square root scale for C >= 0.7
    z <- gamma_rs_sqrt_scale(alpha, gam_param, truncation)
  }
  
  return(z)  # Return the sampled value
}

#### MH sampler for concentration parameter ####
# Function to compute a vectorized sum related to the given parameter alpha
sum_func_vectorized <- function(alpha, n0, log = TRUE) {
  x <- 1:(n0 - 1)  # Create a vector from 1 to n0-1
  a <- x / (alpha + x)^2  # Calculate the function values
  
  max_a <- max(a)  # Find the maximum value of a
  a_prime <- a / max_a  # Normalize a by its maximum value
  
  # Return log of the sum or the sum itself based on the log parameter
  if (log == TRUE) {
    return(log(max_a) + log(sum(a_prime)))  # Logarithmic return
  } else {
    return(sum(x / (alpha + x)^2) )  # Regular sum return
  }
}

# Function to compute the log density of the Jeffrey's prior for alpha
jeffreys_prior_alpha_log_density <- function(a_BNP, pgam) {
  # Return the log density calculation based on the vectorized sum function
  return(0.5 * (sum_func_vectorized(a_BNP, pgam) - log(a_BNP)))
}

# Function to sample new values for alpha based on the Jeffrey's prior
jeffreys_prior_alpha_sampler <- function(a_BNP, tau2_BNP, pgam, k) {
  a_bnp_current <- a_BNP  # Store the current value of a_BNP
  log_a_bnp_current <- log(a_BNP)  # Log of the current value
  
  # Propose a new log value based on a normal distribution
  log_a_bnp_new <- rnorm(1, log_a_bnp_current, sqrt(tau2_BNP))
  a_bnp_new <- exp(log_a_bnp_new)  # Exponentiate to get the new alpha value
  
  # Compute the log posterior for the current and new values
  logpost_a_bnp_current <- lgamma(a_bnp_current) - 
    lgamma(a_bnp_current + pgam) + 
    k * log_a_bnp_current + jeffreys_prior_alpha_log_density(a_bnp_current, pgam)
  
  logpost_a_bnp_new <- lgamma(a_bnp_new) - 
    lgamma(a_bnp_new + pgam) + 
    k * log_a_bnp_new + jeffreys_prior_alpha_log_density(a_bnp_new, pgam)
  
  u <- runif(1)  # Generate a uniform random number for acceptance check
  
  # Calculate the log acceptance probability
  log_accept_prob <- logpost_a_bnp_new - logpost_a_bnp_current + 
    log_a_bnp_new - log_a_bnp_current
  
  # Accept or reject the new value based on the acceptance probability
  if (log(u) <= log_accept_prob) {
    return(a_bnp_new)  # Accept the new value
  } else {
    return(a_bnp_current)  # Reject and return the current value
  }
}


#### main function ####
Blockg.lm <- function(x,y,
                      grp_idx=NULL,
                      adaptive=TRUE,
                      K=40, # Max number of mixture components 
                      burn=10000,
                      nmc=9000,
                      a_BNP=0, # hyperparameter for the nonparametric process (DP vs PY)
                      a_BNP_prior="Gamma", # can be Gamma or Jeffreys
                      # a_BNP==0 implies fully Bayesian with a gamma(1,1) prior on a_BNP
                      a_a_BNP = 1, #hyperparameter for the nonparametric process (DP vs PY)
                      b_a_BNP = 1, #hyperparameter for the nonparametric process (DP vs PY)
                      tau2_BNP=0.05, # variance of proposal for Jeffreys prior sampler 
                      thinning=10, 
                      model.prior="beta-binomial", # uniform, # decaying-beta-binomial 
                      model.init=NULL,
                      hyper.prior="Inv-gamma", # "hyper-g", "hyper-g-n",
                      # "beta-prime-MG","beta-prime"
                      tau2.fixed=FALSE,
                      sigma2.fixed=FALSE,
                      sigma2=1,
                      tau2=1,
                      tau.sampler="slice",
                      hyper.param=NULL,
                      DP.inference=NULL, # can be SB, Dir, Polya-urn
                      prior.corr=TRUE,
                      m=5, # Number of additional values to try to add in Polya urn
                      model=NULL){ 
  
  
  # Initialize the number of observations and predictors
  n <- length(y)            # Number of observations in response variable y
  p <- ncol(x)              # Number of predictors (columns) in matrix x
  
  #### Check for prior correlation and model validity ####
  if (prior.corr == TRUE & !is.null(model)) {
    # Ensure the model size is less than n-1 for valid X'X inverse
    if (sum(model) >= n - 1) {
      stop("X'X inverse in g-prior is not defined for p > n case. 
          Select a model such that model size is < n-1.")
    }
  }
  
  #### Set default inference method for adaptive analysis if not specified ####
  if (adaptive == TRUE & is.null(DP.inference)) {
    DP.inference <- "Polya-urn"  
  }
  
  #### Handle grouping index (grp_idx) ####
  if (is.null(grp_idx)) {
    if (adaptive == FALSE) {
      stop("Either the argument grp_idx should be provided, or adaptive should be set to TRUE")
    }
    
    # Initialize grp_idx based on provided model or default settings
    if (is.null(model)) {
      if (!is.null(model.init)) {
        grp_idx <- model.init
      } else {
        grp_idx <- if (DP.inference == "Polya-urn") rep(0, p) else rep(1, p)
      }
    } else {
      grp_idx <- model  # Use the provided model as grp_idx
    }
    
  } else {
    # Ensure that grp_idx has valid values
    if (adaptive == FALSE) {
      K <- K0 <- length(unique(grp_idx))
    }
    
    # Check the length of grp_idx
    if (length(grp_idx) != p) {
      stop("The argument grp_idx must have length p, where p is the number of columns in X.")
    }
    
    # Validate that grp_idx is an ordered sequence
    if (any(diff(unique(grp_idx)) > 1)) {
      stop("Groups skip a number somewhere in grp_idx. Ensure that grp_idx is an ordered sequence
         from 1 to p with no skips. A valid example is 1,1,1,2,2,3,3,3,4,5,5.")
    }
  }
  
  #### Initialize the random flag for adaptive analysis ####
  random <- FALSE
  if (adaptive & a_BNP == 0) {
    if (a_BNP_prior == "Jeffreys") {
      # Handle case for Jeffreys prior
      if (DP.inference == "SB") {
        stop("Jeffreys prior for a_BNP is currently not implemented with SB inference procedure")
      }
      a_BNP <- rf(1, 1, 1)  # Draw from F(1,1) distribution
    } else if (a_BNP_prior == "Gamma") {
      a_BNP <- rgamma(1, a_a_BNP, b_a_BNP)
    } else {
      stop("a_BNP priors currently allowed are Jeffreys prior or Gamma(a_a_BNP, b_a_BNP) prior")
    }
    random <- TRUE
  }
  
  #### Hyperparameter checks based on specified prior ####
  if (hyper.prior == "Inv-gamma") {
    if (is.null(hyper.param)) {
      a1 <- 1/2
      a2 <- n / 2
    }
    if (length(hyper.param) == 2) {
      if (hyper.param[1] > 0 & hyper.param[2] > 0) {
        a1 <- hyper.param[1]
        a2 <- hyper.param[2]
      } else {
        stop("Inv-gamma requires both parameters > 0")
      }
    }
    if (length(hyper.param) > 2 | length(hyper.param) == 1) {
      stop("Inv-gamma specification requires only two hyper-parameters, both > 0")
    }
  }
  
  if (hyper.prior == "hyper-g") {
    if (is.null(hyper.param)) {
      a <- 3
    } else {
      if (length(hyper.param) > 1 | (hyper.param != 3 & hyper.param != 4)) {
        stop("Hyper-g requires just one hyper parameter with recommended value 3 or 4")
      }
    }
  }
  
  if (hyper.prior == "hyper-g-n") {
    if (is.null(hyper.param)) {
      a <- 3
    } else {
      if (length(hyper.param) > 1 | (hyper.param != 3 & hyper.param != 4)) {
        stop("Hyper-g/n requires just one hyper parameter with value 3 or 4")
      }
    }
  }
  
  if (hyper.prior == "beta-prime-MG") {
    if (is.null(hyper.param)) {
      a <- -3 / 4
    } else {
      if (length(hyper.param) > 1) {
        stop("Maryuma and George recommend a to be between -1 and -1/2 and b to be adaptive")
      }
      if (!(-1 < hyper.param & hyper.param < -1/2)) {
        stop("Maryuma and George recommend a to be between -1 and -1/2")
      }
    }
  }
  
  if (hyper.prior == "beta-prime") {
    if (is.null(hyper.param)) {
      a <- -0.5
      b_param <- -0.5
    }
    if (length(hyper.param) == 2) {
      if (hyper.param[1] > -1 & hyper.param[2] > -1) {
        a <- hyper.param[1]
        b_param <- hyper.param[2]
      } else {
        stop("Beta-prime specification requires a and b > -1")
      }
    }
    if (length(hyper.param) > 2 | length(hyper.param) == 1) {
      stop("Beta-prime specification requires only two hyper-parameters, both > -1")
    }
  }
  
  #### Initialization ####
  # Create matrices to save variables
  GammaSave <- matrix(NA, nmc, p)
  BetaSave <- matrix(NA, nmc, p + 1)
  Sigma2Save <- matrix(NA, nmc, 1)
  Tau2Save <- matrix(NA, nmc, 1)
  logBF212Save <- matrix(NA, nmc, 1)
  g1g2equal2Save <- matrix(NA, nmc, 1)
  nclusterSave <- matrix(NA, nmc, 1)
  aBNPSave <- matrix(NA, nmc, 1)
  aBNPacceptSave <- matrix(NA, nmc, 1)
  logmargSave <- matrix(NA, nmc, 1)
  gvalSave <- matrix(NA, nmc, p)
  grpidSave <- matrix(NA, nmc, p)
  timemat <- matrix(NA, nmc * thinning + burn, 7)
  modelSave <- matrix(NA, nmc, 5)
  
  # Calculate X'X matrix
  xtx.full <- t(x) %*% x
  # Initialize parameters
  g.old <- n  # Previous group size
  
  # Set model parameters
  if (!is.null(model)) {
    if (!is.null(model.init) & !all(model.init == model)) {
      stop("model.init should be the same as model or null when model is specified.")
    }
    gam <- model
    pgam <- sum(gam)  # Sum of selected variables
    b <- rep(0, pgam)  # Initialize coefficients for selected variables
  } else {
    if (!is.null(model.init)) {
      gam <- model.init  # Use initial model if provided
    } else {
      gam <- rep(0, p)  # All variables are excluded by default
    }
    b <- rep(0, p)  # Initialize coefficients for all variables
  }
  
  
  
  # Initialize parameters for the model
  alpha <- mean(y)
  sigma2 <- sigma2
  tau2 <- tau2
  nu <- rinvgamma(1, 1/2, 1)
  
  # Handle clustering parameters
  eta <- rbeta(1, 1, 1)
  grp_idx.gam <- grp_idx[as.logical(gam)]
  group_ids.gam <- unique(grp_idx.gam)
  K0 <- length(group_ids.gam)  # Number of unique groups
  if (adaptive == TRUE) {
    if (DP.inference == "Polya-urn") {
      K <- K0  # Set number of clusters based on unique groups
    }
  }
  
  # Count subjects in each group
  n_k <- count_subjects(K, grp_idx.gam)
  
  # CAUTION
  g_K <- rep(g.old,K) # if adaptive=T and DP.inf="PU" and grp_idx is specified, 
  #technically we should give K unique g_K values but ignoring it for now since
  # we start from null model, so it doesn't matter 
  
  # Adaptive inference for clustering
  if (adaptive == TRUE) {
    if (DP.inference == "SB") {
      stick <- stick_break(K, n_k + 1, a_BNP + sum(n_k) - cumsum(n_k), log = TRUE)
      lclust_prob <- stick$p
    } else if (DP.inference == "Dir") {
      lclust_prob <- log(rgamma(K, a_BNP / K + n_k, 1))
    }
  }
  
  
  
  #### MCMC Iteration loop ####
  for(t in 1:(nmc*thinning+burn)){
    if (t%% 2000 == 0){
      #cat(paste(t," "))
      print(paste("=======Iteration no.", t,"=======\n"))
      print(paste("a_BNP= ", a_BNP,"\n"))
    }
    
    #### Update Gamma####
    start_time <- Sys.time()
    # Check if the model is NULL (i.e., not previously defined)
    if(is.null(model)){
      
      # Propose a new gamma value using the proposal function
      gam.prop <- proposal.gamma(gam)
      
      if(adaptive==TRUE){
        if(DP.inference=="Polya-urn"){
          # Identify differences in the proposed gamma
          gam.diff <- which(gam.prop!=gam)
          # Store current values of gam, grp_idx, g_K etc with suffix .p (initialization for proposed)
          grp_idx.p <- grp_idx
          g_K.p <- g_K
          grp_idx.gam.p <- grp_idx.p[as.logical(gam)]
          
          # Order group indices for proposed gamma
          grp_idx.gam.p <- get_lexi_order2(grp_idx.gam.p)
          group_ids.gam.p = unique(grp_idx.gam.p)
          K.p=K0.p = length(group_ids.gam.p)
          n_k.p <- count_subjects(K.p, grp_idx.gam.p)
          
          # Iterate through differences in gamma
          for(i in 1:length(gam.diff)){
            ind <- gam.diff[i] # Index of the gamma difference
            
            # Case where gamma changes from 0 to 1
            if(gam[ind]==0 & gam.prop[ind]==1){
              
              # If there are no current clusters, start a new cluster
              if(length(n_k.p)==0){
                grp_idx.p[ind] <- 1 # Initialize group index
                
                # sample g_k from prior distribution
                if(hyper.prior=="Inv-gamma"){
                  g_K.new <- rinvgamma(1,a1,a2) 
                }else if(hyper.prior=="hyper-g"){
                  g_K.new <- rhyperg(1,a)
                }else if(hyper.prior=="hyper-g-n"){
                  g_K.new <- rhypergn(1,a,n)
                }else if(hyper.prior=="beta-prime"|hyper.prior=="beta-prime-MG"){
                  if(hyper.prior=="beta-prime"){
                    b_param <- b_param
                  }else if(hyper.prior=="beta-prime-MG"){
                    b_param <- (n-pgam-5)/2 -a
                  }
                  bp.tran <- rbeta(1,a+1,b_param+1)
                  g_K.new <- exp(log(bp.tran)-log(1-bp.tran))
                }
                g_K.p <- g_K.new
                
              }else{
                # If there are existing clusters, Sample from the categorical distribution for group index for the newly added variable 
                grp_idx.p[ind] <- rcatlp(1,log(c(n_k.p,a_BNP)))+1
                
                # If the new group index exceeds current maximum, sample g_K from prior
                if(grp_idx.p[ind]==K.p+1){
                  # sample g_k from prior distribution
                  if(hyper.prior=="Inv-gamma"){
                    g_K.new <- rinvgamma(1,a1,a2) 
                  }else if(hyper.prior=="hyper-g"){
                    g_K.new <- rhyperg(1,a)
                  }else if(hyper.prior=="hyper-g-n"){
                    g_K.new <- rhypergn(1,a,n)
                  }else if(hyper.prior=="beta-prime"|hyper.prior=="beta-prime-MG"){
                    if(hyper.prior=="beta-prime"){
                      b_param <- b_param
                    }else if(hyper.prior=="beta-prime-MG"){
                      b_param <- (n-pgam-5)/2 -a
                    }
                    bp.tran <- rbeta(1,a+1,b_param+1)
                    g_K.new <- exp(log(bp.tran)-log(1-bp.tran))
                  }
                  g_K.p <- c(g_K.p,g_K.new)
                }
              }
              
              # Update group indices based on the proposed gamma
              temp <- grp_idx.p[as.logical(gam.prop)]
              grp_idx.gam.p <- temp[temp>0]
              g_K.p <- g_K.p[unique(grp_idx.gam.p)]
              grp_idx.gam.p <- get_lexi_order2(grp_idx.gam.p)
              temp[temp>0] <- grp_idx.gam.p
              
              # Recalculate unique group IDs and count subjects
              group_ids.gam.p = unique(grp_idx.gam.p)
              K.p=K0.p = length(group_ids.gam.p)
              n_k.p <- count_subjects(K.p, grp_idx.gam.p)
              
              # Update proposed group index
              grp_idx.p <- rep(0,p)
              grp_idx.p[as.logical(gam.prop)] <- temp
              
              
            }else if(gam[ind]==1 & gam.prop[ind]==0){
              # Case where gamma changes from 1 to 0
              temp <- grp_idx.p[as.logical(gam.prop)]
              grp_idx.gam.p <- temp[temp>0]
              
              # Update group indices and counts
              g_K.p <- g_K.p[unique(grp_idx.gam.p)]
              grp_idx.gam.p <- get_lexi_order2(grp_idx.gam.p)
              temp[temp>0] <- grp_idx.gam.p
              
              # Recalculate unique group IDs and count subjects
              group_ids.gam.p = unique(grp_idx.gam.p)
              K.p=K0.p = length(group_ids.gam.p)
              n_k.p <- count_subjects(K.p, grp_idx.gam.p)
              
              # Update proposed group index
              grp_idx.p <- rep(0,p)
              grp_idx.p[as.logical(gam.prop)] <- temp#grp_idx.gam.p
              
              
            }
            
          }
          
        }else{
          # Handle case when DP inference is not Polya-urn
          g_K.p <- g_K
          grp_idx.p <- grp_idx
        }
      }else{
        # If not adaptive, retain previous values
        g_K.p <- g_K
        grp_idx.p <- grp_idx
      }
      
      # Calculate log marginal likelihood for proposed and current models
      if(sum(gam.prop)< n-1){
        logmarg.prop.obj <- log_marginal(x,y,g_K.p,grp_idx.p,gam.prop,hyper.prior=hyper.prior,
                                         prior.corr=prior.corr,tau2)
      }
      logmarg.curr.obj <- log_marginal(x,y,g_K,grp_idx, gam,hyper.prior=hyper.prior,
                                       prior.corr=prior.corr,tau2)

      # Calculate acceptance probability for proposed gamma and delta
      # under beta-binomial
      if(model.prior=="beta-binomial"){
        if(sum(gam.prop)<n-1){
          oj.num.log <- -lchoose(p,sum(gam.prop)) + logmarg.prop.obj$logmg  
        }else{
          oj.num.log <- -Inf
        }
        oj.den.log <-  -lchoose(p,sum(gam))+logmarg.curr.obj$logmg
        
        # logbf21 <- logmarg.m2.obj$logmg - (logmarg.m1.obj$logmg)
      }else if (model.prior=="Uniform"){ # Under Uniform model prior
        # define oj.num.log and oj.den.log
        if(sum(gam.prop)<n-1){
          oj.num.log <- logmarg.prop.obj$logmg  
        }else{
          oj.num.log <- -Inf
        }
        
        oj.den.log <- logmarg.curr.obj$logmg
        # logbf21 <- logmarg.m2.obj$logmg - logmarg.m1.obj$logmg
      }else if(model.prior=="decaying-beta-binomial"){
        a.bb <- 1
        #s0 <- min(n,p)
        #k.bb <- (p-s0/2)/(log(p)*s0/2)
        if(n<=p){b.bb=1}else{
          b.bb <- log(p)#k.bb*log(p)  
        }
        
        # define oj.num.log and oj.den.log
        if(sum(gam.prop)<n-1){
          oj.num.log <- logmarg.prop.obj$logmg + lgamma(a.bb+sum(gam.prop))+lgamma(b.bb+p-sum(gam.prop))#dbbinom(sum(gam.prop), p, a.bb, b.bb,log = TRUE)  
        }else{
          oj.num.log <- -Inf
        }
        
        oj.den.log <- logmarg.curr.obj$logmg + lgamma(a.bb+sum(gam))+lgamma(b.bb+p-sum(gam))#dbbinom(sum(gam), p, a.bb, b.bb,log = TRUE) 
        
      }
      
      if(DP.inference=="Polya-urn"){
        # Compute log density ratio for Jeffreys prior
        if(a_BNP_prior=="Jeffreys"){
          if((sum(gam.prop)==2 & sum(gam)==1)||
             (sum(gam.prop)==1 & sum(gam)==2)||
             (sum(gam.prop)==0 & sum(gam)==1)||
             (sum(gam.prop)==1 & sum(gam)==0)){
            log_a_BNP_density_ratio <- 0
          }else{
            log_a_BNP_density_ratio <- jeffreys_prior_alpha_log_density(a_BNP, sum(gam.prop)) -
              jeffreys_prior_alpha_log_density(a_BNP, sum(gam))
          }
        }else{
          log_a_BNP_density_ratio <-0
        }
      }
      
      # Generate a random number for acceptance check
      u.gam <- runif(1)
      
      if(DP.inference=="Polya-urn"){
        # Calculate log acceptance probability
        loga.gam.prob <- min(as.numeric(oj.num.log-oj.den.log)+log_a_BNP_density_ratio,0)
      }else{
        # Calculate log acceptance probability
        loga.gam.prob <- min(as.numeric(oj.num.log-oj.den.log),0)
        
      }
      
      # Check if the proposed gamma is accepted
      if(log(u.gam)<=loga.gam.prob){
        gam <- gam.prop
        xtx <- logmarg.prop.obj$xtx
        xty <- logmarg.prop.obj$xty
        grp_idx <- logmarg.prop.obj$grp_idx
        g_K <- logmarg.prop.obj$g_K
        ggam <- logmarg.prop.obj$ggam
        sig2rate <- logmarg.prop.obj$sig2rate
        logmarg <- logmarg.prop.obj$logmg
      }else{
        gam <- gam
        xtx <- logmarg.curr.obj$xtx
        xty <- logmarg.curr.obj$xty
        grp_idx <- logmarg.curr.obj$grp_idx
        g_K <- logmarg.curr.obj$g_K
        ggam <- logmarg.curr.obj$ggam
        sig2rate <- logmarg.curr.obj$sig2rate
        logmarg <- logmarg.curr.obj$logmg
      }
      
    }else{
      logmarg.obj <- log_marginal(x,y,g_K,grp_idx,gam,hyper.prior=hyper.prior,
                                  prior.corr=prior.corr,tau2)
      gam <- gam
      xtx <- logmarg.obj$xtx
      xty <- logmarg.obj$xty
      grp_idx <- logmarg.obj$grp_idx
      g_K <- logmarg.obj$g_K
      ggam <- logmarg.obj$ggam
      sig2rate <- logmarg.obj$sig2rate
      logmarg <- logmarg.obj$logmg
      logbf21 <- NA
    }
    
    
    end_time <- Sys.time()
    timemat[t,1] <- end_time-start_time
    
    
    #### Update alpha and Beta #####
    # Start tracking the execution time for the alpha update
    start_time <- Sys.time()
    
    # Sample alpha from a normal distribution with mean as the mean of y and standard deviation adjusted by sigma2 and sample size n
    alpha <- rnorm(1, mean(y), sqrt(sigma2/n))
    
    # Calculate the number of non-zero elements in gam
    pgam <- sum(gam)
    
    # Check if xtx (cross-product matrix) is NULL
    if (is.null(xtx)) {
      # If xtx is NULL, initialize b (regression coefficients) as a zero vector
      b <- rep(0, p)
    } else {
      # Create a diagonal matrix with the inverse square roots of ggam
      ggam_inv_half <- diag(1/sqrt(ggam), nrow = pgam, ncol = pgam)
      
      # Create the precision matrix based on whether prior.corr is TRUE or FALSE
      if (prior.corr == TRUE) {
        # If prior.corr is TRUE, adjust the precision matrix with tau2
        prec.mat <- (ggam_inv_half %*% xtx %*% ggam_inv_half) / tau2 + xtx
      } else {
        # If prior.corr is FALSE, adjust the precision matrix differently
        prec.mat <- (diag(1/ggam, nrow = pgam, ncol = pgam)) / tau2 + xtx
      }
      
      # Compute the posterior distribution for b
      # Covariance matrix is the inverse of (1/sigma2 * prec.mat)
      # Mean is given by the product of the covariance matrix and xty/sigma^2
      b <- custom4(1, 1/sigma2 * xty, 1/sigma2 * prec.mat)
      
    }
    
    # Stop tracking the execution time for alpha update
    end_time <- Sys.time()
    
    # Record the time taken for the alpha update
    timemat[t, 2] <- end_time - start_time
    
    #### Update sigma2 ####
    # Start tracking the execution time for the sigma2 update
    start_time <- Sys.time()
    
    # Check if sigma2 is not fixed
    if (sigma2.fixed == FALSE) {
      # Sample sigma2 from an inverse gamma distribution
      sigma2 <- rinvgamma(1, (n - 1) / 2, sig2rate / 2)
    }
    
    # Stop tracking the execution time for sigma2 update
    end_time <- Sys.time()
    
    # Record the time taken for the sigma2 update
    timemat[t, 3] <- end_time - start_time
    
    
    #### Update grp_idx (1,..p) ####
    # Start tracking the execution time for the adaptive update
    start_time <- Sys.time()
    
    # Check if adaptive inference is enabled
    if (adaptive == TRUE) {
      
      # Check if Polya-urn inference is being used
      if (DP.inference == "Polya-urn") {
        # Identify active groups based on gam
        grp_idx.gam <- grp_idx[as.logical(gam)]
        
        # Initialize cluster counts
        K = K0 = K0_in_mod = length(unique(grp_idx.gam))
        var_in_mod <- which(gam == 1)
        pgam <- length(var_in_mod)  # Number of active variables
        
        # Proceed if there are multiple active variables
        if (pgam > 1) {
          for (r in 1:pgam) {
            ind.var <- var_in_mod[r]  # Current variable index
            grp_idx.gam.n <- grp_idx.gam  # Copy of group indices
            
            # Temporarily remove the current variable from the group index
            grp_idx.gam.n[r] <- 0
            uni_clust_left <- unique(grp_idx.gam[-r])  # Unique clusters left
            g_K.n <- g_K[uni_clust_left]  # Updated cluster parameters
            
            # Order the group index
            grp_idx.gam.n[-r] <- get_lexi_order2(grp_idx.gam[-r])
            K.n <- length(unique(grp_idx.gam.n[-r]))  # Count of clusters
            
            # Initialize the q vector
            q <- rep(0, K.n + m)
            
            # Calculate posterior probabilities if there are clusters left
            if (K.n != 0) {
              q[1:K.n] = sapply(1:K.n, cluster_prop_post, 
                                b, xtx, sigma2, g_K.n, 
                                grp_idx.gam.n, r, prior.corr, tau2)
            }
            
            # Sample new hyperparameters based on the specified prior
            if (hyper.prior == "Inv-gamma") {
              g_K.new <- rinvgamma(m, a1, a2)
            } else if (hyper.prior == "hyper-g") {
              g_K.new <- rhyperg(m, a)
            } else if (hyper.prior == "hyper-g-n") {
              g_K.new <- rhypergn(m, a, n)
            } else if (hyper.prior == "beta-prime" | hyper.prior == "beta-prime-MG") {
              if (hyper.prior == "beta-prime") {
                b_param <- b_param
              } else if (hyper.prior == "beta-prime-MG") {
                b_param <- (n - pgam - 5) / 2 - a
              }
              bp.tran <- rbeta(m, a + 1, b_param + 1)
              g_K.new <- exp(log(bp.tran) - log(1 - bp.tran))
            }
            
            # Combine the current and new hyperparameters
            g_K.temp <- c(g_K.n, g_K.new)
            
            # Update the q vector with likelihoods for new clusters
            q[(K.n + 1):(K.n + m)] <- log(a_BNP / m) + 
              sapply((K.n + 1):(K.n + m), 
                     cluster_prop_llikelihood_simplified, 
                     b, xtx, sigma2, g_K.temp, 
                     grp_idx.gam.n, r, prior.corr, tau2)
            
            # Sample new cluster label
            new_label <- rcatlp(1, log_prob = q) + 1
            gstar.new <- g_K.temp[new_label]
            
            # Handle case where a new cluster is created
            if (new_label > K.n) {
              new_label = K.n + 1
              g_K.n <- c(g_K.n, gstar.new)
            }
            
            # Update group indices and cluster parameters
            grp_idx.gam.n[r] <- new_label
            g_K.n <- g_K.n[unique(grp_idx.gam.n)]
            grp_idx.gam.n <- get_lexi_order2(grp_idx.gam.n)
            
            grp_idx.gam <- grp_idx.gam.n
            g_K <- g_K.n
            K = K0 = length(unique(grp_idx.gam))
          }
          
          # Update grp_idx for inactive variables
          grp_idx <- rep(0, p)
          grp_idx[as.logical(gam)] <- grp_idx.gam
          K0_in_mod <- K
        }
        
      } else {  # Alternative inference method
        # Update grp_idx for inactive variables based on prior probabilities
        grp_idx[gam == 0] <- rcatlp(sum(gam == 0), log_prob = lclust_prob) + 1
        
        g <- g_K[grp_idx]  # Updated cluster parameters
        var_in_mod <- which(gam == 1)
        pgam <- length(var_in_mod)  # Number of active variables
        
        # Proceed if there are active variables
        if (pgam != 0) {
          for (r in 1:pgam) {
            ind.var <- var_in_mod[r]  # Current variable index
            likl_cluster <- sapply(1:K, cluster_prop_llikelihood, 
                                   b, xtx, g, gam, sigma2, g_K, r, prior.corr, tau2)
            
            # Update posterior cluster probabilities
            post_lclust_prob <- lclust_prob + likl_cluster
            grp_idx[ind.var] <- rcatlp(1, log_prob = post_lclust_prob) + 1  # Update cluster index
            g <- g_K[grp_idx]  # Update cluster parameters
          }
        }
        
        # Count unique group indices and update K values
        group_ids = unique(grp_idx)
        K0 = length(group_ids)
        K0_in_mod = length(unique(grp_idx[as.logical(gam)]))
        
        # Update clust_prob based on stick-breaking process or gamma representation
        n_k <- count_subjects(K, grp_idx[as.logical(gam)])
        
        if (DP.inference == "SB") {
          stick <- stick_break(K, n_k + 1, a_BNP + sum(n_k) - cumsum(n_k), log = TRUE)
          lclust_prob <- stick$p
        } else {
          lclust_prob <- log(rgamma(K, a_BNP / K + n_k, 1))
        }
      }
      
      # Update a_BNP if random sampling is enabled
      if (random == TRUE) {
        if (DP.inference == "SB") {
          a_BNP <- sample_Plambda_posterior(1, stick$q, alpha0 = a_a_BNP, beta0 = b_a_BNP)
        } else if (DP.inference == "Dir" | DP.inference == "Polya-urn") {
          if (a_BNP_prior == "Gamma") {
            if (pgam != 0) {
              # Sample eta for updating a_BNP
              eta <- rbeta(1, a_BNP + 1, pgam)
              z1 <- rcatlp(1, log_prob = c(log(a_a_BNP + K0_in_mod - 1), 
                                           log(pgam) + log(b_a_BNP - log(eta))))
              if (z1 == 0) {
                a_BNP <- rgamma(1, a_a_BNP + K0_in_mod, b_a_BNP - log(eta))
              } else {
                a_BNP <- rgamma(1, a_a_BNP + K0_in_mod - 1, b_a_BNP - log(eta))
              }
            } else {
              a_BNP <- rgamma(1, a_a_BNP, b_a_BNP)
            }
            aBNPaccept <- 1
          } else if (a_BNP_prior == "Jeffreys") {
            if (pgam > 1) {
              a_BNP_cur <- a_BNP
              a_BNP <- jeffreys_prior_alpha_sampler(a_BNP, tau2_BNP, pgam, K0_in_mod)
              aBNPaccept <- ifelse(a_BNP_cur != a_BNP, 1, 0)
            } else {
              a_BNP <- rf(1, 1, 1)  # Sample from Jeffrey's prior
              aBNPaccept <- 1
            }
          }
        }
      }
    }
    
    # Stop tracking the execution time for the adaptive update
    end_time <- Sys.time()
    
    # Record the time taken for the adaptive update
    timemat[t, 4] <- end_time - start_time
    
    
    
    #### Update g ####
    start_time <- Sys.time()
    # Calculate the number of the active variables
    pgam <- sum(gam)
    
    # Check if the DP inference method is not set to "Polya-urn"
    if(is.null(DP.inference) || DP.inference !="Polya-urn"){
      
      # Extract group-specific parameters
      g <- g_K[grp_idx]
      # Initialize the covariance matrix C based on the input xtx
      if(is.null(xtx)){
        C <- NULL
      }else{
        # Create a vector of coefficients
        bvec <- as.vector(b)
        bvec.full <- rep(0,p)
        bvec.full[which(gam==1)] <- bvec
        # Compute the covariance matrix C depending on whether prior correlation is used
        if(prior.corr==TRUE){
          C <- (2*sigma2*tau2)^{-1}* diag(bvec.full,nrow=length(bvec.full),ncol = 
                                            length(bvec.full)) %*% xtx.full %*% 
            diag(bvec.full, nrow=length(bvec.full),ncol = length(bvec.full))
          
        }else{
          C <- (2*sigma2*tau2)^{-1}* diag(bvec.full^2,nrow=length(bvec.full),ncol = 
                                            length(bvec.full))
        }
      }
      
      # Loop over each group
      for(k in 1:K){
        # Get members of current group
        grp_mem <- which(grp_idx==k)
        grp_mem_in_mod <- which(grp_idx==k & gam==1)
        s_k <- length(grp_mem_in_mod) # size of current group 
        # Check if there are any members in the current group that are also selected
        if(s_k!=0){
          p_k <- length(grp_mem_in_mod)

          # Calculate c_k and d_k based on the current group's covariance
          if(pgam==p_k){
            c_k <- sum(C)
            d_k <- 0
          }else{
            c_k <- sum(C[grp_mem_in_mod,grp_mem_in_mod])
            d_k <- sum(2*C[grp_mem_in_mod,-grp_mem_in_mod]%*% 
                         diag(1/sqrt(g[-grp_mem_in_mod]),nrow = 
                                p-length(grp_mem_in_mod)))
          }
          
          # sample g_k from that distribution using slice + rejection sampling
          if(hyper.prior=="Inv-gamma"){
            # draw from an extended Gamma distribution
            gam_prime <- d_k/(2*sqrt(c_k+a2))
            gk_inv <- ext_gamma_sampler(alpha = a1+p_k/2, gam_param = 
                                          gam_prime,truncation = NULL)
            gk_inv <- gk_inv/(c_k+a2)
            g_K[k] <- 1/gk_inv
            g[grp_mem] <- 1/gk_inv
          }else if(hyper.prior=="hyper-g"){
            t_i <- 1/g_K[k]
            # Slice sampling approach
            ht_i <- (1+t_i)^{-a/2}
            u_i <- runif(1,0,ht_i)
            truncation <- c_k*(u_i^{-2/a}-1)
            zeta <- d_k/(2*sqrt(c_k))
            z_i <- ext_gamma_sampler(alpha = (a+p_k-2)/2,gam_param = zeta,
                                     truncation = truncation)
            gk_inv <- z_i/c_k
            g_K[k] <- 1/gk_inv
            g[grp_mem] <- 1/gk_inv
            
          }else if(hyper.prior=="hyper-g-n"){
            t_i <- 1/g_K[k] 
            # Slice sampling approach
            ht_i <- (1+n*t_i)^{-a/2}
            u_i <- runif(1,0,ht_i)
            truncation <- c_k*1/n*(u_i^{-2/a}-1)
            zeta <- d_k/(2*sqrt(c_k))
            z_i <- ext_gamma_sampler(alpha = (a+p_k-2)/2,gam_param = zeta,
                                     truncation = truncation)
            gk_inv <- z_i/c_k
            g_K[k] <- 1/gk_inv
            g[grp_mem] <- 1/gk_inv
          }else if(hyper.prior=="beta-prime"|hyper.prior=="beta-prime-MG"){
            t_i <- 1/g_K[k] 
            
            # Slice sampling approach
            pgam <- sum(gam)
            if(hyper.prior=="beta-prime-MG"){
              b_param <- (n-pgam-5)/2 -a
            }else if(hyper.prior=="beta-prime"){
              b_param <- b_param
            }
            if(c_k==0 & d_k==0){
              bp.tran <- rbeta(1,a+p_k/2+1,b_param+1-p_k/2)
              g[ind] <- exp(log(bp.tran)-log(1-bp.tran))
            }else{
              ht_i <- (1+t_i)^{-(a+b_param+2)}#{-(n-pgam-1)/2}
              u_i <- runif(1,0,ht_i)
              truncation <- c_k*(u_i^{-1/(a+b_param+2)}-1)#(u_i^{-2/(n-pgam-1)}-1)
              zeta <- d_k/(2*sqrt(c_k))
              z_i <- ext_gamma_sampler(alpha = a+p_k/2+1,gam_param = zeta,
                                       truncation = truncation)
              gk_inv <- z_i/c_k
              g_K[k] <- 1/gk_inv
              g[grp_mem] <- 1/gk_inv
            }
          }}else{
            # Sample g_k from prior distribution if no members in the current group are selected            if(hyper.prior=="Inv-gamma"){
              g_K[k] <- rinvgamma(1,a1,a2) 
            }else if(hyper.prior=="hyper-g"){
              g_K[k] <- rhyperg(1,a)
            }else if(hyper.prior=="hyper-g-n"){
              g_K[k] <- rhypergn(1,a,n)
            }else if(hyper.prior=="beta-prime"|hyper.prior=="beta-prime-MG"){
              if(hyper.prior=="beta-prime"){
                b_param <- b_param
              }else if(hyper.prior=="beta-prime-MG"){
                b_param <- (n-pgam-5)/2 -a
              }
              bp.tran <- rbeta(1,a+1,b_param+1)
              g_K[k] <- exp(log(bp.tran)-log(1-bp.tran))#rbetapr(1, a+1, b_param+1)
            }
            # Update group members with the sampled value
            if(length(grp_mem)!=0){
              g[grp_mem] <- g_K[k]
            }
          }
      
    }else{
      # Handle case where DP inference is set to "Polya-urn"
      if(is.null(xtx)){
        C <- NULL
      }else{
        bvec <- as.vector(b)
        if(prior.corr==TRUE){
          C <- (2*sigma2*tau2)^{-1}* diag(bvec,nrow=length(bvec),ncol = 
                                            length(bvec)) %*% xtx %*% 
            diag(bvec, nrow=length(bvec),ncol = length(bvec))
        }else{
          C <- (2*sigma2*tau2)^{-1}* diag(bvec^2,nrow=length(bvec),ncol = 
                                            length(bvec))
        }
      }
      
      if(K!=0){
        # Loop over each group for Polya-urn
        for(k in 1:K){
          
          g <- g_K[grp_idx.gam]
          
          grp_mem_in_mod <- which(grp_idx.gam==k)
          p_k <- length(grp_mem_in_mod)
          # Calculate c_k and d_k
          if(pgam==p_k){
            c_k <- sum(C)
            d_k <- 0
          }else{
            c_k <- sum(C[grp_mem_in_mod,grp_mem_in_mod])
            d_k <- sum(2*C[grp_mem_in_mod,-grp_mem_in_mod]%*% 
                         diag(1/sqrt(g[-grp_mem_in_mod]),nrow = 
                                pgam-length(grp_mem_in_mod)))
          }
          
          # sample g_k from that distribution using slice + rejection sampling
          if(hyper.prior=="Inv-gamma"){
            # draw from an extended Gamma distribution
            gam_prime <- d_k/(2*sqrt(c_k+a2))
            gk_inv <- ext_gamma_sampler(alpha = a1+p_k/2, gam_param = 
                                          gam_prime,truncation = NULL)
            gk_inv <- gk_inv/(c_k+a2)
            g_K[k] <- 1/gk_inv
          }else if(hyper.prior=="hyper-g"){
            t_i <- 1/g_K[k] 
            ht_i <- (1+t_i)^{-a/2}
            u_i <- runif(1,0,ht_i)
            truncation <- c_k*(u_i^{-2/a}-1)
            zeta <- d_k/(2*sqrt(c_k))
            z_i <- ext_gamma_sampler(alpha = (a+p_k-2)/2,gam_param = zeta,
                                     truncation = truncation)
            gk_inv <- z_i/c_k
            g_K[k] <- 1/gk_inv

          }else if(hyper.prior=="hyper-g-n"){
            t_i <- 1/g_K[k] 
            # Slice sampling approach
            ht_i <- (1+n*t_i)^{-a/2}
            u_i <- runif(1,0,ht_i)
            truncation <- c_k*1/n*(u_i^{-2/a}-1)
            zeta <- d_k/(2*sqrt(c_k))
            z_i <- ext_gamma_sampler(alpha = (a+p_k-2)/2,gam_param = zeta,
                                     truncation = truncation)
            gk_inv <- z_i/c_k
            g_K[k] <- 1/gk_inv
          }else if(hyper.prior=="beta-prime"|hyper.prior=="beta-prime-MG"){
            t_i <- 1/g_K[k] 
            
            # Slice sampling approach
            pgam <- sum(gam)
            if(hyper.prior=="beta-prime-MG"){
              b_param <- (n-pgam-5)/2 -a
            }else if(hyper.prior=="beta-prime"){
              b_param <- b_param
            }
            if(c_k==0 & d_k==0){
              bp.tran <- rbeta(1,a+p_k/2+1,b_param+1-p_k/2)
              g_K[k] <- exp(log(bp.tran)-log(1-bp.tran))
            }else{
              ht_i <- (1+t_i)^{-(a+b_param+2)}#{-(n-pgam-1)/2}
              u_i <- runif(1,0,ht_i)
              truncation <- c_k*(u_i^{-1/(a+b_param+2)}-1)#(u_i^{-2/(n-pgam-1)}-1)
              zeta <- d_k/(2*sqrt(c_k))
              z_i <- ext_gamma_sampler(alpha = a+p_k/2+1,gam_param = zeta,
                                       truncation = truncation)
              gk_inv <- z_i/c_k
              g_K[k] <- 1/gk_inv
              #g[grp_mem] <- 1/gk_inv
            }
          }
        }
      }
      
    }
    
    end_time <- Sys.time()
    timemat[t,5] <- end_time-start_time
    
    #### Update tau2 and nu ####
    if(tau2.fixed==FALSE){
      start_time <- Sys.time()
      if(tau.sampler=="augment"){
        if(pgam>0){
          if(DP.inference=="Polya-urn"){
            ggam <- g_K[grp_idx.gam]
          }else{
            ggam <- g[as.logical(gam)]
            
          }
          ggam_inv_half <- diag(1/sqrt(ggam),nrow = pgam,ncol = pgam)
          if(prior.corr==TRUE){
            tau2rate <- (t(b)%*%(ggam_inv_half %*% xtx %*% ggam_inv_half)%*% b)/(2*sigma2)
          }else{
            tau2rate <- (t(b)%*%(diag(1/ggam,nrow = pgam,ncol = pgam))%*% b)/(2*sigma2)
          }
          #print(tau2)
          tau2 <- 1/rgamma(1, (pgam+1)/2, 1/nu+ tau2rate)
        }else{
          tau2 <- 1/rgamma(1, 1/2,1/nu)
        }
      }else if(tau.sampler=="slice"){
        #### Alternate slice sampler 
        if(pgam>0){
          
          if(is.null(DP.inference) || DP.inference !="Polya-urn"){
            ggam <- g[as.logical(gam)]
          }else{
            ggam <- g_K[grp_idx.gam]
          }
          
          ggam_inv_half <- diag(1/sqrt(ggam),nrow = pgam,ncol = pgam)
          if(prior.corr==TRUE){
            tau2rate <- (t(b)%*%(ggam_inv_half %*% xtx %*% ggam_inv_half)%*% b)/(2*sigma2)
          }else{
            tau2rate <- (t(b)%*%(diag(1/ggam,nrow = pgam,ncol = pgam))%*% b)/(2*sigma2)
          }
          #print(tau2)
          tau2inv <- 1/tau2
          ut <- runif(1, 0,1/(1+tau2inv))
          tau2inv <- rtgamma(1, shape = (pgam+1)/2,rate =  tau2rate, truncation =  1/ut-1)
          tau2 <- 1/tau2inv
          
          #tau2 <- 1/rgamma(1, (pgam+1)/2, 1/nu+ tau2rate)
        }else{
          nu <- 1/rgamma(1, 1/2, 1)
          tau2 <- 1/rgamma(1, 1/2,1/nu)
        }
      }
      
      
      end_time <- Sys.time()
      timemat[t,6] <- end_time-start_time
      
      start_time <- Sys.time()
      
      if(tau.sampler=="augment"){
        nu <- 1/rgamma(1, 1, 1+1/tau2)
      }
      
      end_time <- Sys.time()
      timemat[t,7] <- end_time-start_time
      
    }
    
    #### Case when p=2 ####
    if(p==2){
      m2 <- c(1,1)
      m1 <- c(1,0)
      
      grp_idx.gam <- grp_idx[as.logical(gam)]
      
      
      if(length(grp_idx.gam)==2){
        g_K.temp <- g_K
        grp_idx.temp <- grp_idx
      }else if(length(grp_idx.gam)==1){
        
        new_grpid <- rcatlp(1,log(c(1,a_BNP)))+1
        grp_idx.temp <- grp_idx
        grp_idx.temp[!as.logical(gam)] <- new_grpid
        g_K.temp <- g_K
        
        if(new_grpid==2){
          # sample g_k from prior distribution
          if(hyper.prior=="Inv-gamma"){
            g_K.new <- rinvgamma(1,a1,a2) 
          }else if(hyper.prior=="hyper-g"){
            g_K.new <- rhyperg(1,a)
          }else if(hyper.prior=="hyper-g-n"){
            g_K.new <- rhypergn(1,a,n)
          }else if(hyper.prior=="beta-prime"|hyper.prior=="beta-prime-MG"){
            if(hyper.prior=="beta-prime"){
              b_param <- b_param
            }else if(hyper.prior=="beta-prime-MG"){
              b_param <- (n-pgam-5)/2 -a
            }
            bp.tran <- rbeta(1,a+1,b_param+1)
            g_K.new <- exp(log(bp.tran)-log(1-bp.tran))#rbetapr(1, a+1, b_param+1)
          }
          g_K.temp <- c(g_K.temp,g_K.new)
          
          g_K.temp <- g_K.temp[unique(grp_idx.temp)]
          grp_idx.temp <- get_lexi_order2(grp_idx.temp)
        }
      }else{
        grp_idx.temp <- rep(0,2)
        K.temp <- length(unique(grp_idx.temp[as.logical(gam)]))
        n_k.temp <- count_subjects(K.temp,grp_idx.temp[as.logical(gam)])
        for(h in 1:2){
          
          if(length(n_k.temp)==0){
            grp_idx.temp[h] <- 1
            # sample g_k from prior distribution
            if(hyper.prior=="Inv-gamma"){
              g_K.new <- rinvgamma(1,a1,a2) 
            }else if(hyper.prior=="hyper-g"){
              g_K.new <- rhyperg(1,a)
            }else if(hyper.prior=="hyper-g-n"){
              g_K.new <- rhypergn(1,a,n)
            }else if(hyper.prior=="beta-prime"|hyper.prior=="beta-prime-MG"){
              if(hyper.prior=="beta-prime"){
                b_param <- b_param
              }else if(hyper.prior=="beta-prime-MG"){
                b_param <- (n-pgam-5)/2 -a
              }
              bp.tran <- rbeta(1,a+1,b_param+1)
              g_K.new <- exp(log(bp.tran)-log(1-bp.tran))#rbetapr(1, a+1, b_param+1)
            }
            g_K.temp <- g_K.new
            n_k.temp <- 1
            K.temp <- 1
            
          }else{
            
            grp_idx.temp[h] <- rcatlp(1,log(c(n_k.temp,a_BNP)))+1
            
            if(grp_idx.temp[h]==K.temp+1){
              # sample g_k from prior distribution
              if(hyper.prior=="Inv-gamma"){
                g_K.new <- rinvgamma(1,a1,a2) 
              }else if(hyper.prior=="hyper-g"){
                g_K.new <- rhyperg(1,a)
              }else if(hyper.prior=="hyper-g-n"){
                g_K.new <- rhypergn(1,a,n)
              }else if(hyper.prior=="beta-prime"|hyper.prior=="beta-prime-MG"){
                if(hyper.prior=="beta-prime"){
                  b_param <- b_param
                }else if(hyper.prior=="beta-prime-MG"){
                  b_param <- (n-pgam-5)/2 -a
                }
                bp.tran <- rbeta(1,a+1,b_param+1)
                g_K.new <- exp(log(bp.tran)-log(1-bp.tran))#rbetapr(1, a+1, b_param+1)
              }
              g_K.temp <- c(g_K.temp,g_K.new)
            }
          }
        }
      }
      
      logmarg.m2.obj <- log_marginal(x,y,g_K.temp, grp_idx.temp,m2,hyper.prior=hyper.prior,
                                     prior.corr=prior.corr,tau2)
      logmarg.m1.obj <- log_marginal(x,y,g_K.temp,grp_idx.temp,m1,hyper.prior=hyper.prior,
                                     prior.corr=prior.corr,tau2)
      logbf21 <- logmarg.m2.obj$logmg - (logmarg.m1.obj$logmg)
      g1_g2_equal <- as.numeric(grp_idx.temp[1]==grp_idx.temp[2])
      
      temp.vec <- rep(0,5)
      
      if(setequal(gam,c(0,0))){
        temp.vec[1]=1
      }else if(setequal(gam,c(1,0))){
        temp.vec[2]=1
      }else if(setequal(gam,c(0,1))){
        temp.vec[3]=1
      }else if(setequal(gam,c(1,1))){
        if(grp_idx.temp[1]==grp_idx.temp[2]){
          temp.vec[4]=1
        }else{
          temp.vec[5]=1
        }
      }
      
    }else{
      logbf21 <- NA
      g1_g2_equal <- NA
      temp.vec <- rep(NA,5)
    }
    
    
    
    #### Store results as a list ####
    betastore=rep(0,p)
    #Save results post burn-in
    if(t > burn && (t - burn) %% thinning == 0){
      rr  <- floor((t - burn) / thinning)
      GammaSave[rr, ] = gam
      betastore[which(gam==1)] <- b
      betastore <- c(alpha,betastore)
      BetaSave[rr, ] = betastore
      Sigma2Save[rr, ] <- sigma2
      Tau2Save[rr, ] <- tau2
      if(adaptive){
        nclusterSave[rr, ] <- K0_in_mod
      }else{
        nclusterSave[rr, ] <- K
      }
      logmargSave[rr, ] <- logmarg
      aBNPSave[rr, ] <- a_BNP
      if(random){
        aBNPacceptSave[rr, ] <- aBNPaccept
      }else{
        aBNPacceptSave[rr, ] <- NA
      }
      if(is.null(DP.inference) || DP.inference !="Polya-urn"){
        gvalSave[rr, ] <- g
      }else{
        gvalSave[rr, as.logical(gam)] <- g_K[grp_idx]
      }
      grpidSave[rr, ] <- grp_idx
      logBF212Save[rr, ] <- logbf21
      g1g2equal2Save[rr, ] <- g1_g2_equal
      modelSave[rr, ] <- temp.vec
    }
  }
  
  
  result <- list("BetaSamples"=BetaSave,
                 "GammaSamples"=GammaSave,
                 "Sigma2Samples"=Sigma2Save,
                 "Tau2Samples"=Tau2Save,
                 "ncluster"=nclusterSave,
                 "aBNPSamples"=aBNPSave,
                 "aBNPaccept"=aBNPacceptSave,
                 "grpid"=grpidSave,
                 "gsamples"=gvalSave,
                 "logBF21"=logBF212Save,
                 "gequal"= g1g2equal2Save,
                 "modelfreq"=modelSave,
                 "timemat"=timemat)
  # Return object
  return(result)
  
}



