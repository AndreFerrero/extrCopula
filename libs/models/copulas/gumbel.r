source("libs/packages.R")

copula_gumbel <- list(

  name = "gumbel",

  simulate_u = function(theta, n) {

  # --- Independence case ---
  if (theta == 1) {
    return(runif(n))
  }

  # --- Invalid parameter ---
  if (theta < 1) return(NULL)

  # --- Gumbel Archimedean simulation ---
  val_gamma <- (cos(pi / (2 * theta)))^theta

  V <- stabledist::rstable(
    n = 1,
    alpha = 1 / theta,
    beta  = 1,
    gamma = val_gamma,
    delta = 0,
    pm    = 1
  )

  if (!is.finite(V) || V <= 0) return(NULL)

  E <- rexp(n)
  exp(-(E / V)^(1 / theta))
  },

  log_density = function(u, theta) {
    copula::dCopula(u, copula::gumbelCopula(theta, dim = length(u)), log = TRUE)
  },

  log_prior = function(theta, a = 2, b = 1) {
    if (theta <= 1) return(-Inf)
    dgamma(theta - 1, a, b, log = TRUE)
  },

  diag = function(u, theta, n = length(u)) {
    
    if (theta < 1) return(NA)             # invalid parameter
    if (theta == 1) return(u^n)           # independence case
    
    u^(n^(1/theta))
}
)
