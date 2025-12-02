library(copula)

set.seed(123)

# -----------------------------------------------------------
# 1. TRUE parameter and dimension
# -----------------------------------------------------------
n <- 100
theta_true <- 3

# -----------------------------------------------------------
# 2. Simulate ONE 100-dimensional vector from Gumbel copula
# -----------------------------------------------------------
cop <- gumbelCopula(param = theta_true, dim = n)
U <- as.numeric(rCopula(1, cop)) # true uniforms

# -----------------------------------------------------------
# 3. Generate X via lognormal margins
# -----------------------------------------------------------
mu <- 0
sigma <- 1
X <- qlnorm(U, meanlog = mu, sdlog = sigma)

# -----------------------------------------------------------
# 4A. Parametric margins
# -----------------------------------------------------------
mu_hat <- mean(log(X))
sigma_hat <- sd(log(X))
u_hat_param <- plnorm(X, meanlog = mu_hat, sdlog = sigma_hat)

# -----------------------------------------------------------
# 4B. ECDF pseudo-observations
# -----------------------------------------------------------
u_hat_ecdf <- rank(X) / (n + 1)

# -----------------------------------------------------------
# 5. Gumbel log-likelihood
# -----------------------------------------------------------
loglik_gumbel <- function(theta, u) {
  if (theta <= 1) {
    return(-1e10)
  } # keep optimizer stable
  cop <- gumbelCopula(param = theta, dim = length(u))
  ll <- log(dCopula(u, copula = cop))
  return(ll)
}

# Negative log-likelihood for optimization
negLL <- function(theta, u) -loglik_gumbel(theta, u)

# -----------------------------------------------------------
# 6. Compute estimators via optim()
# -----------------------------------------------------------

# TRUE UNIFORMS
est_trueU <- optim(
  par = 5, fn = negLL, u = U, method = "L-BFGS-B",
  lower = 1.001, upper = 30
)$par

# PARAMETRIC MARGINS
est_param <- optim(
  par = 5, fn = negLL, u = u_hat_param, method = "L-BFGS-B",
  lower = 1.001, upper = 30
)$par

# ECDF PSEUDO-OBSERVATIONS
est_ecdf <- optim(
  par = 5, fn = negLL, u = u_hat_ecdf, method = "L-BFGS-B",
  lower = 1.001, upper = 30
)$par

# -----------------------------------------------------------
# 7. Likelihood curves for comparison
# -----------------------------------------------------------
theta_grid <- seq(1.01, 12, length.out = 200)

ll_trueU <- sapply(theta_grid, loglik_gumbel, u = U)
ll_param <- sapply(theta_grid, loglik_gumbel, u = u_hat_param)
ll_ecdf <- sapply(theta_grid, loglik_gumbel, u = u_hat_ecdf)

plot(theta_grid, ll_trueU,
  type = "l", lwd = 3, col = "black",
  xlab = expression(theta), ylab = "Log-likelihood",
  main = "Gumbel log-likelihood shapes"
)

lines(theta_grid, ll_param, col = "blue", lwd = 2)
lines(theta_grid, ll_ecdf, col = "red", lwd = 2)

abline(v = theta_true, col = "darkgreen", lty = 2, lwd = 2)
abline(v = est_trueU, col = "black", lty = 3, lwd = 2)
abline(v = est_param, col = "blue", lty = 3, lwd = 2)
abline(v = est_ecdf, col = "red", lty = 3, lwd = 2)

legend("topright",
  legend = c(
    "True U", "Parametric Margin", "ECDF",
    "True theta", "theta_trueU", "theta_param", "theta_ecdf"
  ),
  col = c(
    "black", "blue", "red", "darkgreen",
    "black", "blue", "red"
  ),
  lty = c(1, 1, 1, 2, 3, 3, 3),
  lwd = c(3, 2, 2, 2, 2, 2, 2)
)

# -----------------------------------------------------------
# 8. Print estimates
# -----------------------------------------------------------
cat("\nESTIMATED THETAS:\n")
cat("theta_true =", theta_true, "\n")
cat("theta_trueU   =", est_trueU, "\n")
cat("theta_param   =", est_param, "\n")
cat("theta_ecdf    =", est_ecdf, "\n")
