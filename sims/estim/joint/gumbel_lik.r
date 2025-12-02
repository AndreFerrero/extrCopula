library(copula)

# -----------------------------------------------------------
# 1. Choose dimension and true parameter
# -----------------------------------------------------------
n <- 100            # dimension of the sample vector
alpha_true <- 2    # true Gumbel parameter

# -----------------------------------------------------------
# 2. Simulate ONE vector from the n-dimensional Gumbel copula
# -----------------------------------------------------------
cop <- gumbelCopula(param = alpha_true, dim = n)
U <- as.numeric(rCopula(1, cop))   # one vector of uniforms

# -----------------------------------------------------------
# 3. Choose margins and generate the observed sample X
# -----------------------------------------------------------
# lognormal margins just for concreteness:
mu <- 0
sigma <- 1
X <- qlnorm(U, meanlog = mu, sdlog = sigma)

# -----------------------------------------------------------
# 4. Estimate margins and obtain pseudo-observations
# -----------------------------------------------------------
mu_hat <- mean(log(X))
sigma_hat <- sd(log(X))

u_hat <- plnorm(X, meanlog = mu_hat, sdlog = sigma_hat)

# -----------------------------------------------------------
# 5. Define log-likelihood for the n-dimensional Gumbel copula
# -----------------------------------------------------------
loglik_gumbel <- function(alpha, u) {
  if (alpha <= 1) return(-Inf)
  cop <- gumbelCopula(param = alpha, dim = length(u))
  ll <- log(dCopula(u, copula = cop))
  return(ll)
}

# -----------------------------------------------------------
# 6. Evaluate likelihood over a grid of alpha values
# -----------------------------------------------------------
alpha_grid <- seq(1.0, 5, length.out = 200)
ll_values <- sapply(alpha_grid, loglik_gumbel, u = u_hat)

# -----------------------------------------------------------
# 7. Plot the likelihood curve
# -----------------------------------------------------------
plot(alpha_grid, ll_values, type = "l", lwd = 2,
     xlab = expression(alpha),
     ylab = "Log-likelihood",
     main = paste("Log-likelihood for ONE sample vector (n =", n, ")"))
abline(v = alpha_true, col = "red", lwd = 2, lty = 2)
legend("bottomright", legend = c("True alpha"), col = "red", lwd = 2, lty = 2)
