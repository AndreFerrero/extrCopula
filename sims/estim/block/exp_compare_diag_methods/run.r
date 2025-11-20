library(copula)
plots_folder <- here("sims", "estim", "block", "exp_compare_diag_methods", "results", "plots")

# ----------------------------
# Simulation parameters
# ----------------------------
set.seed(46)
nk <- 12                  # block size (dimension of copula)
theta_true <- 4          # true Gumbel copula parameter
K_vec <- c(5, 10, 20, 50, 100) # number of blocks
nSim <- 100             # number of simulation replications

# ----------------------------
# Functions
# ----------------------------

# Simulate block-dependent data
simulate_blocks <- function(K, cop) {
  lapply(1:K, function(b) as.numeric(rCopula(1, cop)))
}

# Analytical closed-form DMLE
theta_closed <- function(blocks, n) {
  Y <- sapply(blocks, max)
  theta <- log(n) / (log(length(blocks)) - log(sum(-log(Y))))
  return(max(min(theta, 10), 1))
}

# Numerical diagonal MLE
psi_gumbel <- function(t, theta) {
  exp(-t^(1/theta))
}

psi_inv_gumbel <- function(u, theta) {
  (-log(u))^theta
}

psi_prime_gumbel <- function(t, theta) {
  (-1/theta) * t^(1/theta - 1) * exp(-t^(1/theta))
}

psi_inv_prime_gumbel <- function(u, theta) {
  -theta/u * (-log(u))^(theta - 1)
}

diag_prime_gumbel <- function(u, n, theta) {
  inv_val  <- psi_inv_gumbel(u, theta)        # ψ⁻¹(u)
  psi_p    <- psi_prime_gumbel(n * inv_val,   # ψ'(n ψ⁻¹(u))
                               theta)
  inv_p    <- psi_inv_prime_gumbel(u, theta)  # (ψ⁻¹)'(u)

  diag_prime <- log(n * psi_p * inv_p)
  return(diag_prime)
}

diagLogLik_max <- function(theta, blocks, n) {
  ll <- 0
  for(block in blocks) {
    ll <- ll + diag_prime_gumbel(max(block), n, theta)
  }
  return(-ll)
}

theta_numerical <- function(blocks, n) {
  res <- optim(par=1.5, fn=diagLogLik_max, blocks=blocks, n=n,
               method="L-BFGS-B", lower=1, upper=10)
  res$par
}

# ----------------------------
# Likelihood shape
# ----------------------------

theta_grid <- seq(1.1, 10, length.out = 200)

# Choose sets of simulated blocks
blocks_test_theta1p5 <- simulate_blocks(K = 1000, cop = gumbelCopula(1.5, dim = nk))
blocks_test_theta4   <- simulate_blocks(K = 1000, cop = gumbelCopula(4,   dim = nk))
blocks_test_theta6   <- simulate_blocks(K = 1000, cop = gumbelCopula(6,   dim = nk))

# Compute log-likelihood curves
dloglik_theta1p5 <- sapply(theta_grid, function(th) -diagLogLik_max(th, blocks_test_theta1p5, nk))
dloglik_theta4   <- sapply(theta_grid, function(th) -diagLogLik_max(th, blocks_test_theta4,   nk))
dloglik_theta6   <- sapply(theta_grid, function(th) -diagLogLik_max(th, blocks_test_theta6,   nk))

png(file = here(plots_folder, "dloglik_theta.png"), width = 1200, height = 400)
par(mfrow = c(1,3))

# θ = 1.5
plot(theta_grid, dloglik_theta1p5, type="l", lwd=2,
     xlab=expression(theta), ylab="Log-likelihood", main="Theta = 1.5")
abline(v = 1.5, col="red", lty=2)

# θ = 4
plot(theta_grid, dloglik_theta4, type="l", lwd=2,
     xlab=expression(theta), ylab="Log-likelihood", main="Theta = 4")
abline(v = 4, col="red", lty=2)

# θ = 6
plot(theta_grid, dloglik_theta6, type="l", lwd=2,
     xlab=expression(theta), ylab="Log-likelihood", main="Theta = 6")
abline(v = 6, col="red", lty=2)

par(mfrow = c(1,1))
dev.off()


# ----------------------------
# Simulation study
# ----------------------------
results_closed <- matrix(NA, nrow=nSim, ncol=length(K_vec))
results_numerical <- matrix(NA, nrow=nSim, ncol=length(K_vec))
colnames(results_closed) <- colnames(results_numerical) <- paste0("K=", K_vec)

cop <- gumbelCopula(theta_true, dim=nk)

for(j in seq_along(K_vec)) {
  K <- K_vec[j]
  for(i in 1:nSim) {
    blocks <- simulate_blocks(K, nk, cop)
    results_closed[i,j] <- theta_closed(blocks, nk)
    results_numerical[i,j] <- theta_numerical(blocks, nk)
  }
  cat("Completed simulations for K =", K, "\n")
}

# ----------------------------
# Summary statistics
# ----------------------------
# means_closed <- colMeans(results_closed)
# means_numerical <- colMeans(results_numerical)
# sds_closed <- apply(results_closed, 2, sd)
# sds_numerical <- apply(results_numerical, 2, sd)

# print(data.frame(K=K_vec,
#                  mean_closed=means_closed, sd_closed=sds_closed,
#                  mean_numerical=means_numerical, sd_numerical=sds_numerical))

# ----------------------------
# Visualization
# ----------------------------
png(file = here(plots_folder, "compare_gumbel_theta4.png"))
par(mfrow=c(1,2))
boxplot(results_closed, names=paste0("K=", K_vec),
        main="Closed-form DMLE", col="lightblue", ylab=expression(hat(theta)))
abline(h=theta_true, col="red", lwd=2, lty=2)

boxplot(results_numerical, names=paste0("K=", K_vec),
        main="Numerical DMLE", col="lightgreen", ylab=expression(hat(theta)))
abline(h=theta_true, col="red", lwd=2, lty=2)
par(mfrow = c(1,1))
dev.off()
