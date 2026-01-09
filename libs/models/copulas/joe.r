# ------------------------------
# Joe Copula Object
# ------------------------------
copula_joe <- list(

  name = "joe",

  # --------------------------
  # Generator and inverse generator
  # --------------------------
  psi = function(t, theta) 1 - (1 - exp(-t))^(1/theta),   # generator
  inv_psi = function(u, theta) -log(1 - (1 - u)^theta),

  # --------------------------
  # 1. Simulate uniforms using latent variable method
  # --------------------------
  simulate_u = function(theta, n) {

    # Independence
    if (theta == 1) return(runif(n))
    if (theta < 1) return(NULL)  # invalid

    # --------------------------
    # 1a. Compute discrete latent PMF
    # P(V = k) = (-1)^(k+1) * choose(1/theta, k)
    # Truncate at K_max
    # --------------------------
    K_max <- 1000
    probs <- sapply(1:K_max, function(k) { (-1)^(k+1) * choose(1/theta, k) })
    probs <- probs[probs > 0]
    probs <- probs / sum(probs)  # normalize to sum 1

    # --------------------------
    # 1b. Sample latent V
    # --------------------------
    V <- sample(1:length(probs), size = 1, prob = probs)

    # --------------------------
    # 1c. Generate uniforms
    # U_i = 1 - (1 - exp(-E_i / V))^1/theta
    # --------------------------
    E <- rexp(n)
    U <- copula_joe$psi(E/V, theta)
    U
  },

  # --------------------------
  # 2. Log-density using copula package
  # --------------------------
  log_density = function(u, theta) {
    copula::dCopula(u, copula::joeCopula(theta, dim = length(u)), log = TRUE)
  },

  # --------------------------
  # 3. Log-prior for theta
  # --------------------------
  log_prior = function(theta, a = 2, b = 1) {
    if (theta <= 1) return(-Inf)
    dgamma(theta - 1, shape = a, rate = b, log = TRUE)
  },

  # --------------------------
  # 4. Diagonal function
  # --------------------------
  diag = function(u, theta, n = length(u)) {

    if (theta < 1) return(NA)   # invalid
    if (theta == 1) return(u^n) # independence

    # Joe copula diagonal in n dimensions: C(u,...,u)
    # Approximation: C(u,...,u) = 1 - (1 - u)^(n^(1/theta))
    1 - (1 - u)^(n^(1/theta))
  }
)
