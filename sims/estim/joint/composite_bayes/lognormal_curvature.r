library(copula)
library(here)
library(coda)
library(mvtnorm)
library(future)
library(future.apply)
library(numDeriv) # Required for Gradient/Hessian
library(psych)

set.seed(123)

# -----------------------------------------------------------
# 1. Simulate Data
# -----------------------------------------------------------
n <- 1000
theta_true <- 3
mu_true <- 0
sigma_true <- 1

cat("1. Simulating Data...\n")
cop <- gumbelCopula(param = theta_true, dim = n)
U <- as.numeric(rCopula(1, cop))
X <- qlnorm(U, meanlog = mu_true, sdlog = sigma_true)

# -----------------------------------------------------------
# 2. Setup Pairwise Indices
# -----------------------------------------------------------
# We use a subset of pairs to keep the "Composite" part dominant but manageable
n_pairs_keep <- 3000
all_pairs <- t(combn(n, 2))
keep_idx <- sample(nrow(all_pairs), n_pairs_keep)
pair_indices <- all_pairs[keep_idx, ]

# -----------------------------------------------------------
# 3. Helper Functions & Transformations
# -----------------------------------------------------------

# We work in Transformed Space to ensure unconstrained optimization
# 1. Theta -> lambda = log(theta - 1)
# 2. Mu    -> mu (No transform needed, defined on Real line)
# 3. Sigma -> zeta = log(sigma)

to_natural <- function(vec) {
  theta <- exp(vec[1]) + 1
  mu    <- vec[2]
  sigma <- exp(vec[3])
  return(c(theta, mu, sigma))
}

to_transformed <- function(vec) {
  lambda <- log(vec[1] - 1)
  mu     <- vec[2]
  zeta   <- log(vec[3])
  return(c(lambda, mu, zeta))
}

# The Unadjusted Composite Log-Likelihood
# (Marginal Likelihood + Pairwise Copula Likelihood)
calc_unadjusted_ll <- function(trans_vec, data, pairs_idx) {
  nat <- to_natural(trans_vec)
  theta <- nat[1]; mu <- nat[2]; sigma <- nat[3]
  
  # Bounds check (for safety during optim)
  if(theta > 20 || sigma <= 0.01) return(-1e10)
  
  # 1. Marginals
  d_marg <- dlnorm(data, meanlog = mu, sdlog = sigma, log = TRUE)
  ll_marg <- sum(d_marg)
  
  # 2. Copula (Pairwise)
  u_hat <- plnorm(data, meanlog = mu, sdlog = sigma)
  u_hat <- pmin(pmax(u_hat, 1e-9), 1 - 1e-9)
  
  u_pairs <- cbind(u_hat[pairs_idx[,1]], u_hat[pairs_idx[,2]])
  d_cop <- dCopula(u_pairs, copula = gumbelCopula(theta, dim = 2))
  ll_cop <- sum(log(pmax(d_cop, 1e-100)))
  
  return(ll_marg + ll_cop)
}

# The Prior Function (Unadjusted)
calc_priors <- function(trans_vec) {
  nat <- to_natural(trans_vec)
  theta <- nat[1]; mu <- nat[2]; sigma <- nat[3]
  
  lp <- 0
  # Theta: Gamma(2,1)
  lp <- lp + dgamma(theta - 1, 2, 1, log = TRUE)
  # Mu: Normal(0, 10)
  lp <- lp + dnorm(mu, 0, 10, log = TRUE)
  # Sigma: LogNormal(0, 1) - Boundary Avoiding
  lp <- lp + dlnorm(sigma, 0, 1, log = TRUE)
  
  # Jacobian Adjustments for the transformations
  # d(theta)/d(lambda) = theta - 1
  lp <- lp + log(theta - 1)
  # d(sigma)/d(zeta) = sigma
  lp <- lp + log(sigma)
  
  return(lp)
}

# -----------------------------------------------------------
# IMPROVED: Weighted composite likelihood
# -----------------------------------------------------------
calc_unadjusted_ll_weighted <- function(trans_vec, data, pairs_idx) {
  nat <- to_natural(trans_vec)
  theta <- nat[1]; mu <- nat[2]; sigma <- nat[3]
  
  if(theta <= 1.001 || theta > 50 || sigma <= 0.001) return(-1e20)
  
  # 1. Marginals (Weight = 1)
  d_marg <- dlnorm(data, meanlog = mu, sdlog = sigma, log = TRUE)
  if(any(!is.finite(d_marg))) return(-1e20)
  ll_marg <- sum(d_marg)
  
  # 2. Copula (Weighted)
  # Calculate Scaling Factor: 
  # We want Pairwise contribution to be roughly same magnitude as Marginals
  n_data <- length(data)
  n_pairs <- nrow(pairs_idx)
  w_pair <- n_data / n_pairs
  
  u_hat <- plnorm(data, meanlog = mu, sdlog = sigma)
  eps <- 1e-6
  u_hat <- pmin(pmax(u_hat, eps), 1 - eps)
  
  u_pairs <- cbind(u_hat[pairs_idx[,1]], u_hat[pairs_idx[,2]])
  
  d_cop <- dCopula(u_pairs, copula = gumbelCopula(theta, dim = 2), log = TRUE)
  d_cop[!is.finite(d_cop)] <- -1e20
  
  ll_cop <- sum(d_cop)
  
  # APPLY WEIGHT
  # This balances the "force" of the gradients
  return(ll_marg + (w_pair * ll_cop))
}

# -----------------------------------------------------------
# SURFACE INSPECTION TOOL
# -----------------------------------------------------------
inspect_surface <- function(data, pairs, center_param, range_mult = 0.5, n_grid = 40) {
  
  # Center param should be in TRANSFORMED scale
  # e.g., c(log(theta-1), mu, log(sigma))
  
  # Define Ranges based on center +/- range_mult
  lam_vec <- seq(center_param[1] - range_mult, center_param[1] + range_mult, length.out = n_grid)
  mu_vec  <- seq(center_param[2] - range_mult, center_param[2] + range_mult, length.out = n_grid)
  zet_vec <- seq(center_param[3] - range_mult, center_param[3] + range_mult, length.out = n_grid)
  
  # -------------------------------------------------------
  # PLOT 1: Lambda (Theta) vs Zeta (Sigma) -- Fixing Mu
  # -------------------------------------------------------
  z_matrix <- matrix(NA, n_grid, n_grid)
  
  cat("Calculating Grid 1: Theta vs Sigma...\n")
  for(i in 1:n_grid) {
    for(j in 1:n_grid) {
      # Fix Mu at center, vary others
      p <- c(lam_vec[i], center_param[2], zet_vec[j]) 
      z_matrix[i,j] <- calc_unadjusted_ll_weighted(p, data, pairs)
    }
  }
  
  # Convert axes back to natural for labeling
  theta_axis <- exp(lam_vec) + 1
  sigma_axis <- exp(zet_vec)
  
  # Identify Max in this slice
  max_idx <- which(z_matrix == max(z_matrix), arr.ind=TRUE)
  
  filled.contour(theta_axis, sigma_axis, z_matrix,
                 main = "LL Surface: Theta vs Sigma (Mu Fixed)",
                 xlab = "Theta", ylab = "Sigma",
                 color.palette = terrain.colors,
                 plot.axes = {
                   axis(1); axis(2);
                   points(exp(center_param[1])+1, exp(center_param[3]), pch=19, col="red", cex=1.5);
                   points(theta_axis[max_idx[1]], sigma_axis[max_idx[2]], pch=4, col="blue", cex=1.5, lwd=2)
                 },
                 key.title = title(main="LogLik"))
  
  cat("Red Dot = Start/Center, Blue X = Grid Maximum\n")
}

# -----------------------------------------------------------
# RUN DIAGNOSTIC
# -----------------------------------------------------------

# 1. Define where you want to look (Truth or Current Mode)
# Let's look around the TRUE values first to see if the global max is where we expect
truth_trans <- to_transformed(c(theta_true, mu_true, sigma_true))

# 2. Run inspection
# Adjust range_mult to zoom in/out (0.5 means +/- 0.5 in log scale)
inspect_surface(X, pair_indices, truth_trans, range_mult = 1.0, n_grid = 50)

# -----------------------------------------------------------
# 4. CALIBRATION (Updated for Weighted Likelihood)
# -----------------------------------------------------------
cat("2. Finding Mode and Calibrating Curvature (Weighted)...\n")

# A. Optimization using the WEIGHTED function
#    (We define the objective function wrapper)
obj_fn_weighted <- function(p) {
  # Negate because optim minimizes
  -1 * (calc_unadjusted_ll_weighted(p, X, pair_indices) + calc_priors(p))
}

# Run Optim (BFGS is best for this smooth surface)
# Start close to truth or use the two-stage approach discussed earlier
start_trans <- to_transformed(c(1.5, mean(log(X)), sd(log(X)))) 
opt <- optim(start_trans, obj_fn_weighted, method = "BFGS", hessian = FALSE)

mode_trans <- opt$par
cat("   Mode found (Natural):", round(to_natural(mode_trans), 4), "\n")

# -----------------------------------------------------------
# B. Calculate H (Hessian of the WEIGHTED posterior)
# -----------------------------------------------------------
# We calculate the Hessian of the objective function at the mode.
# Note: obj_fn_weighted is (-LL), so the Hessian is positive definite.
# We want H = -Hessian(LL), so H_hat = Hessian(obj_fn) is actually correct 
# because H should be the Information matrix (Positive Definite).
# Let's stick to the notation: H = - Hessian(LogPost).

H_hat <- hessian(function(p) calc_unadjusted_ll_weighted(p, X, pair_indices) + calc_priors(p), mode_trans)
H_hat <- -H_hat # Convert to Negative Hessian (should be Neg Definite)
H_hat <- -H_hat # Flip again to get Information Matrix (Pos Definite for computations)
# Wait, simpler: H_hat should be the Observed Information Matrix (Positive Definite)
H_hat <- hessian(obj_fn_weighted, mode_trans) 

# -----------------------------------------------------------
# C. Calculate J (Variance of the WEIGHTED Score)
# -----------------------------------------------------------
# This is the most critical part. We must simulate data and 
# calculate the Gradient of the WEIGHTED likelihood.

n_boot <- 200 # Increased for stability
grad_matrix <- matrix(NA, nrow=n_boot, ncol=3)
mode_nat <- to_natural(mode_trans)

# Generative model at the estimated mode
sim_cop <- gumbelCopula(mode_nat[1], dim = n)

cat("   Bootstrapping for J matrix (Weighted)...\n")
for(i in 1:n_boot) {
  # 1. Simulate Dataset
  U_sim <- as.numeric(rCopula(1, sim_cop))
  X_sim <- qlnorm(U_sim, meanlog = mode_nat[2], sdlog = mode_nat[3])
  
  # 2. Calculate Gradient of the WEIGHTED Likelihood
  #    Must use the SAME weights as optimization
  g <- grad(function(p) calc_unadjusted_ll_weighted(p, X_sim, pair_indices), mode_trans)
  
  grad_matrix[i, ] <- g
}

J_hat <- cov(grad_matrix)

# -----------------------------------------------------------
# D. Compute Adjustment Matrix C
# -----------------------------------------------------------
# Calculate target variance (Godambe)
Cov_target <- solve(H_hat) %*% J_hat %*% solve(H_hat)

# Calculate current implied variance (Inverse Hessian)
Cov_current <- solve(H_hat)

# Compute C using Cholesky decomposition (Chandler & Bate)
M_H <- chol(H_hat)        
M_G <- chol(solve(Cov_target)) 

C_adj <- solve(M_H) %*% M_G 

cat("   Adjustment Matrix C (Diagonal):\n")
print(diag(C_adj))

# -----------------------------------------------------------
# 5. Adjusted Posterior Function for MCMC
# -----------------------------------------------------------

log_posterior_adjusted <- function(param_vec, data, pairs, mode_p, C_mat) {
  # param_vec is the MCMC proposal (in the "Broad" space)
  
  # 1. WARP PARAMETERS toward the mode (shrinkage)
  # theta* = mode + C (theta - mode)
  # If C is small (0.1), a large step in param_vec becomes a small step in theta*
  # This makes the likelihood surface "look" wider to the sampler.
  diff <- param_vec - mode_p
  param_star <- mode_p + as.vector(C_mat %*% diff) # or t(C_mat)? C is usually symmetric-ish
  
  # 2. Evaluate Unadjusted Likelihood at Warped Parameters
  ll_val <- calc_unadjusted_ll(param_star, data, pairs)
  
  # 3. Add Priors (Evaluated at the PROPOSED/BROAD parameters, not warped)
  # Rationale: The prior describes our belief about the parameters *before* data.
  # The adjustment corrects the *Likelihood*.
  lp_val <- calc_priors(param_vec) 
  
  return(ll_val + lp_val)
}

# -----------------------------------------------------------
# 6. Worker Function
# -----------------------------------------------------------
run_worker_adj <- function(data, n_iter, init, chain_id, pairs, mode_p, C_mat) {
  library(copula); library(mvtnorm); library(coda)
  
  chain <- matrix(NA, n_iter, 3)
  colnames(chain) <- c("lambda", "mu", "zeta")
  curr <- init
  
  # Wrapper
  lp_fn <- function(p) log_posterior_adjusted(p, data, pairs, mode_p, C_mat)
  
  curr_lp <- lp_fn(curr)
  if(!is.finite(curr_lp)) stop("Bad Init")
  
  # Adaptive Settings
  # Since we adjusted the curvature, the posterior should look roughly 
  # like the prior or a standard normal covariance structurally.
  # We still use adaptation to fine-tune.
  cov_mat <- diag(0.01, 3); sd_scale <- 2.38^2/3; eps <- diag(1e-6, 3)
  
  for(i in 1:n_iter) {
    prop <- as.numeric(rmvnorm(1, curr, cov_mat))
    prop_lp <- lp_fn(prop)
    
    if(is.finite(prop_lp) && log(runif(1)) < (prop_lp - curr_lp)) {
      curr <- prop; curr_lp <- prop_lp
    }
    chain[i,] <- curr
    
    if(i > 500 && i %% 50 == 0) cov_mat <- sd_scale * cov(chain[1:i,]) + eps
    
    # Progress
    if(i %% 1000 == 0) {
      status_file <- file.path(tempdir(), paste0("status_", chain_id, ".txt"))
      writeLines(as.character(i/n_iter), status_file)
    }
  }
  return(mcmc(chain))
}

# -----------------------------------------------------------
# 7. Execution
# -----------------------------------------------------------
n_chains <- 3
plan(multisession, workers = n_chains)

# Start MCMC near the mode (since we know where it is)
inits_list <- list()
for(c in 1:n_chains) inits_list[[c]] <- mode_trans + rnorm(3, 0, 0.1)

cat("3. Launching Adjusted MCMC...\n")
futures <- list()
for(c in 1:n_chains) {
  futures[[c]] <- future({
    run_worker_adj(X, 10000, inits_list[[c]], c, pair_indices, mode_trans, C_adj)
  }, seed=TRUE)
}

while(!all(resolved(futures))) Sys.sleep(1)
cat("Done.\n")

# -----------------------------------------------------------
# 8. Analysis
# -----------------------------------------------------------
chains_raw <- value(futures)

# Transform chains back to Natural Scale
chains_nat <- lapply(chains_raw, function(ch) {
  # The chain contains the "Broad" (Adjusted) parameters.
  # These ARE the parameters of interest (with correct variance).
  # We just transform them from log-space to natural space.
  t(apply(ch, 1, to_natural))
})

mcmc_clean <- window(mcmc.list(lapply(chains_nat, mcmc)), start=2000, thin=5)
colnames(mcmc_clean[[1]]) <- c("theta", "mu", "sigma")

cat("\n--- DIAGNOSTICS ---\n")
print(effectiveSize(mcmc_clean))

# Visuals
par(mfrow=c(1,3))
plot(mcmc_clean[, "theta"], main="Theta", density=FALSE, auto.layout=FALSE)
plot(mcmc_clean[, "sigma"], main="Sigma", density=FALSE, auto.layout=FALSE)
plot(mcmc_clean[, "mu"], main="Mu", density=FALSE, auto.layout=FALSE)

cat("\nTrue Values: Theta=3, Mu=0, Sigma=1\n")
cat("Posterior Means:\n")
print(colMeans(as.matrix(mcmc_clean)))

# Check correlation recovery
pairs.panels(as.matrix(mcmc_clean))