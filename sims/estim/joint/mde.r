# ===============================================================
# Single-Realization Archimedean Copula Estimation (theta, v, h)
# Supports "clayton" and "gumbel" families
# Uses L2 minimum distance between empirical CDF and model G_v
# ===============================================================

# ---------- Helper functions ----------

# Inverse generator function psi^{-1}(u) for given family and theta
psi_inv <- function(u, family, theta) {
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)  # avoid boundaries
  switch(family,
         clayton = u^(-theta) - 1,
         gumbel  = (-log(u))^theta,
         stop("Unsupported copula family"))
}

# Gaussian kernel CDF estimator for marginal F
kernel_cdf <- function(x_eval, data, h) {
  m <- outer(x_eval, data, FUN = function(x, z) pnorm((x - z)/h))
  rowMeans(m)
}

# Empirical CDF at a grid of points
empirical_cdf <- function(x_eval, data) {
  sapply(x_eval, function(t) mean(data <= t))
}

# ---------- Objective function builder ----------
# Returns a function of par = c(log_theta, log_h, log_v)
make_loss_fn <- function(data, family, x_grid) {
  function(par) {
    theta <- exp(par[1]); h <- exp(par[2]); v <- exp(par[3])
    
    # Smoothed marginal CDF
    Fhat <- kernel_cdf(x_grid, data, h)
    Fhat <- pmin(pmax(Fhat, 1e-12), 1 - 1e-12)
    
    # Model CDF G_v
    Gv <- exp(-v * psi_inv(Fhat, family, theta))
    
    # Empirical CDF
    Fn <- empirical_cdf(x_grid, data)
    
    # L2 distance (Riemann approximation)
    dx <- if(length(x_grid) > 1) diff(x_grid)[1] else 1
    loss <- sum((Fn - Gv)^2) * dx
    
    if(!is.finite(loss)) loss <- 1e12
    loss
  }
}

# ---------- Simulate single vector ----------
simulate_vector <- function(n, family, theta, margin = "norm") {
  if(!requireNamespace("copula", quietly = TRUE)) install.packages("copula")
  library(copula)
  
  cop <- switch(family,
                clayton = claytonCopula(param = theta, dim = n),
                gumbel  = gumbelCopula(param = theta, dim = n),
                stop("Unsupported family"))
  
  U <- as.numeric(rCopula(1, cop))
  
  switch(margin,
         norm = qnorm(U),
         exp  = qexp(U, rate = 1),
         U)  # uniform if unspecified
}

# ---------- Main estimation procedure ----------
estimate_copula <- function(x, family) {
  n <- length(x)
  x_grid <- seq(min(x) - 0.2, max(x) + 0.2, length.out = 200)
  
  loss_fn <- make_loss_fn(x, family, x_grid)
  
  # Initial parameters (log-scale)
  init_par <- c(log(1.0), log(0.5 * sd(x)), log(1.0))
  
  # Lower/upper bounds (log-scale)
  lower <- c(log(1e-3), log(1e-4), log(1e-6))
  upper <- c(log(1e2),  log(5 * sd(x) + 1e-6), log(1e6))
  
  # Optimize using nlminb
  opt <- nlminb(start = init_par, objective = loss_fn,
                lower = lower, upper = upper,
                control = list(eval.max = 5000, iter.max = 5000, rel.tol = 1e-8))
  
  est <- exp(opt$par)
  names(est) <- c("theta", "h", "v")
  
  # Fitted CDF for diagnostics
  Fhat <- kernel_cdf(x_grid, x, est["h"])
  Ghat <- exp(-est["v"] * psi_inv(Fhat, family, est["theta"]))
  Fn <- empirical_cdf(x_grid, x)
  
  list(opt = opt, est = est, x_grid = x_grid, Fn = Fn, Ghat = Ghat)
}

# ---------- Example / demonstration ----------
set.seed(42)
n <- 600
true_family <- "gumbel"
true_theta <- 2.0

# Simulate data
x <- simulate_vector(n, true_family, true_theta, margin = "norm")

# Estimate copula parameters
res <- estimate_copula(x, true_family)
cat(sprintf("Estimated theta = %.4f (true %.4f)\nEstimated h = %.4f\nEstimated v = %.4f\n",
            res$est["theta"], true_theta, res$est["h"], res$est["v"]))

# ---------- Diagnostic plot ----------
plot(res$x_grid, res$Fn, lwd = 2, ylim = c(0,1),
     main = sprintf("Empirical CDF vs Fitted G_v (family=%s)", true_family),
     ylab = "CDF", xlab = "x")
lines(res$x_grid, res$Ghat, lwd = 2, col = "red")
legend("bottomright", legend = c("Empirical CDF", "Fitted G_v"),
       col = c("black", "red"), lwd = 2)
