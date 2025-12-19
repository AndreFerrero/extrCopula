source("libs/packages.R")

# ------------------------------
# Clayton copula functions
# ------------------------------

copula_clayton <- list(

  name = "clayton",

  # --------------------------
  # 1. Simulate uniforms using latent gamma variable
  # --------------------------
  # Clayton copula CDF: C(u1,...,un) = (sum(u_i^-theta - 1) + 1)^(-1/theta)
  # Latent variable V ~ Gamma(1/theta, 1)
  # Then U_i = (1 + E_i / V)^(-1/theta), E_i ~ Exp(1)
  simulate_u = function(theta, n) {

  # --- Independence case ---
  if (abs(theta) < 1e-8) {
    return(runif(n))
  }

  # --- Invalid parameter ---
  if (theta < 0) return(NULL)

  # --- Clayton Archimedean simulation ---
  V <- rgamma(n = 1, shape = 1 / theta, rate = 1)
  if (!is.finite(V) || V <= 0) return(NULL)

  E <- rexp(n)
  (1 + E / V)^(-1 / theta)
},

  # --------------------------
  # 2. Log-density of the Clayton copula
  # --------------------------
  log_density = function(u, theta) {
    copula::dCopula(u, copula::claytonCopula(theta, dim = length(u)), log = TRUE)
  },

  # --------------------------
  # 3. Log-prior for the copula parameter
  # --------------------------
  log_prior = function(theta, a = 2, b = 1) {
    if (theta <= 0) return(-Inf)
    dgamma(theta, shape = a, rate = b, log = TRUE)
  },

  diag = function(u, theta) {
    n <- length(u)
    
    if (theta < 0) return(NA)             # invalid parameter
    if (theta == 0) return(u^n)           # independence case
    
    (n * (u^(-theta) - 1) + 1)^(-1/theta)
  }
) 
