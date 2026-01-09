# ------------------------------
# t Copula Object (Fast Version)
# ------------------------------
copula_t <- list(
  
  name = "t",
  
  # --------------------------
  # 1. Simulate uniforms along the diagonal
  # --------------------------
  simulate_u = function(theta, n, df = 4) {
    # theta: equicorrelation along the diagonal (0 ≤ theta ≤ 1)
    # df: degrees of freedom of t-copula
    if(theta < 0 || theta > 1) stop("theta must be in [0,1]")
    
    # Equicorrelation matrix
    Sigma <- matrix(theta, n, n)
    diag(Sigma) <- 1
    
    # Cholesky decomposition for MVN
    L <- chol(Sigma)
    
    # Z ~ MVN(0, Sigma)
    Z <- L %*% rnorm(n)
    
    # W ~ Chi-squared(df)
    W <- rchisq(1, df)
    
    # Multivariate t vector
    X <- Z / sqrt(W / df)
    
    # Transform to uniforms via t CDF
    U <- pt(X, df)
    as.numeric(U)
  },
  
  # --------------------------
  # 2. Log-density along the diagonal
  # --------------------------
  log_density = function(u, theta, df = 4) {
    # Exact log density (diagonal vector)
    n <- length(u)
    Sigma <- matrix(theta, n, n)
    diag(Sigma) <- 1
    
    # Transform uniforms to t quantiles
    x <- qt(u, df)
    
    # Log-density of multivariate t
    mvtnorm::dmvt(x, sigma = Sigma, df = df, log = TRUE) - sum(dt(x, df, log = TRUE))
  },
  
  # --------------------------
  # 3. Prior for copula parameter
  # --------------------------
  log_prior = function(theta, a = 2, b = 2) {
    # Beta prior on [0,1]
    if(theta < 0 || theta > 1) return(-Inf)
    dbeta(theta, a, b, log = TRUE)
  },
  
  # --------------------------
  # 4. Diagonal function for maxima
  # --------------------------
  diag = function(u, theta, n = length(u)) {
    # Approximate "power distortion" along the diagonal for maxima
    # Crude formula: effective correlation scales the power
    # For independent case (theta=0) u^n
    # For perfect dependence (theta=1) u
    u^((1 - theta) * n + theta)
  }
)