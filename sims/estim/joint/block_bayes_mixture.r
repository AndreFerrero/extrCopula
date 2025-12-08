library(copula)
library(coda)
library(mvtnorm) # Needed for multivariate normal proposals
library(future) 
library(future.apply)

set.seed(123)

# -----------------------------------------------------------
# 1. Simulate Data
# -----------------------------------------------------------
n <- 150 # Borderline for full likelihood, but manageable
theta_true <- 3
mu_true <- 0
sigma_true <- 1

cat("1. Simulating Data (n =", n, ")...\n")
cop <- gumbelCopula(param = theta_true, dim = n)
U <- as.numeric(rCopula(1, cop))
X <- qlnorm(U, meanlog = mu_true, sdlog = sigma_true)
Y <- log(X)

# -----------------------------------------------------------
# 2. Model Definitions
# -----------------------------------------------------------
d_mix <- function(y, w, m1, s1, m2, s2) {
  w * dnorm(y, m1, s1) + (1 - w) * dnorm(y, m2, s2)
}

p_mix <- function(y, w, m1, s1, m2, s2) {
  w * pnorm(y, m1, s1) + (1 - w) * pnorm(y, m2, s2)
}

log_posterior <- function(param_vec, data_y) {
  # Unpack
  theta <- param_vec[1]; w <- param_vec[2]
  mu1 <- param_vec[3]; s1 <- param_vec[4]
  mu2 <- param_vec[5]; s2 <- param_vec[6]

  # --- PRIORS ---
  if (theta <= 1.01 || theta > 20) return(-Inf)
  # Gamma prior on theta (shifted)
  lp_theta <- dgamma(theta - 1, shape = 2, rate = 1, log = TRUE)

  if (w <= 0.01 || w >= 0.99) return(-Inf)
  if (mu1 > mu2) return(-Inf) # Label switching constraint
  if (s1 <= 0.05 || s2 <= 0.05) return(-Inf) # Singularity constraint

  # Priors on marginals
  lp_s <- dgamma(s1, 3, 2, log = TRUE) + dgamma(s2, 3, 2, log = TRUE)
  lp_mu <- dnorm(mu1, 0, 5, log = TRUE) + dnorm(mu2, 0, 5, log = TRUE)

  # --- MARGINAL LIKELIHOOD ---
  dens_vals <- d_mix(data_y, w, mu1, s1, mu2, s2)
  if (any(dens_vals <= 1e-300)) return(-Inf)
  ll_margin <- sum(log(dens_vals))

  # --- COPULA LIKELIHOOD ---
  u_hat <- p_mix(data_y, w, mu1, s1, mu2, s2)
  # Clamp for numerical safety
  u_hat <- pmin(pmax(u_hat, 1e-7), 1 - 1e-7)

  # FULL Likelihood (n-dimensional)
  # Note: dCopula returns log density directly with log=TRUE.
  # Valid log-densities are negative, so we ONLY check for -Inf or NaN.
  ld_cop <- dCopula(u_hat, copula = gumbelCopula(theta, dim = length(data_y)), log = TRUE)
  
  if (is.na(ld_cop) || ld_cop == -Inf) return(-Inf)

  return(lp_theta + lp_s + lp_mu + ll_margin + ld_cop)
}

# -----------------------------------------------------------
# 3. Improved Worker: Block Adaptive Metropolis
# -----------------------------------------------------------
run_block_adaptive_chain <- function(data, n_iter, init_vals, chain_id, progress_dir) {
  # Libraries must be loaded inside worker
  library(copula); library(mvtnorm); library(coda)
  
  p_names <- c("theta", "w", "mu1", "s1", "mu2", "s2")
  n_par <- length(p_names)
  
  chain <- matrix(NA, nrow = n_iter, ncol = n_par)
  colnames(chain) <- p_names

  curr_par <- init_vals
  
  # Robust Start: Ensure we don't start at -Inf
  curr_lp <- log_posterior(curr_par, data)
  if (!is.finite(curr_lp)) {
    # Try jittering to find a valid spot
    for(k in 1:20) {
      curr_par <- init_vals * runif(n_par, 0.9, 1.1)
      if(curr_par[3] > curr_par[5]) { # Fix ordering if jitter broke it
         tmp <- curr_par[3]; curr_par[3] <- curr_par[5]; curr_par[5] <- tmp
      }
      curr_lp <- log_posterior(curr_par, data)
      if(is.finite(curr_lp)) break
    }
    if(!is.finite(curr_lp)) stop("Chain failed to initialize: Log-Posterior is -Inf.")
  }

  # --- ADAPTIVE SETTINGS ---
  # Initial small covariance
  cov_mat <- diag(rep(0.005, n_par)) 
  
  # Optimal scaling factor for d=6 (Haario et al)
  sd_scale <- (2.38^2) / n_par 
  eps <- 1e-6 * diag(n_par) # Regularization
  
  adapt_start <- 500 # When to start learning covariance
  total_accepts <- 0
  status_file <- file.path(progress_dir, paste0("status_", chain_id, ".txt"))

  for (i in 1:n_iter) {
    
    # 1. BLOCK PROPOSAL (All params at once)
    # Propose from Multivariate Normal based on current covariance
    prop_val <- rmvnorm(1, mean = curr_par, sigma = cov_mat)[1, ]
    
    prop_lp <- log_posterior(prop_val, data)
    ratio <- prop_lp - curr_lp

    # 2. METROPOLIS STEP
    if (is.finite(ratio) && log(runif(1)) < ratio) {
      curr_par <- prop_val
      curr_lp <- prop_lp
      total_accepts <- total_accepts + 1
    }
    chain[i, ] <- curr_par

    # 3. ADAPTATION (Learn the "Banana Shape")
    if (i > adapt_start && i %% 50 == 0) {
      # Compute empirical covariance of history
      emp_cov <- cov(chain[1:i, ])
      # Update proposal covariance
      cov_mat <- (sd_scale * emp_cov) + eps
    }

    # 4. DASHBOARD
    if (i %% 50 == 0 || i == n_iter) {
      rate <- (total_accepts / i) * 100
      pct <- as.integer((i / n_iter) * 100)
      try({ writeLines(sprintf("%d|%.1f", pct, rate), status_file) }, silent = TRUE)
    }
  }
  return(mcmc(chain))
}

# -----------------------------------------------------------
# 4. Execution
# -----------------------------------------------------------
n_chains <- 3
plan(multisession, workers = n_chains) # Match workers to chains

cat("2. Launching Adaptive Block MCMC...\n")

# Use slightly longer run for Adaptive method to settle
n_iter <- 15000 
burn_in <- 3000

# Prepare Initial Values
inits_list <- list()
for (c in 1:n_chains) {
  inits_list[[c]] <- c(
    theta = runif(1, 2, 4), 
    w = 0.5,
    mu1 = mean(Y) - 0.2, s1 = sd(Y),
    mu2 = mean(Y) + 0.2, s2 = sd(Y)
  )
}

# Setup Status Files
progress_dir <- tempdir()
file.remove(list.files(progress_dir, pattern = "status_", full.names = TRUE))

# Launch Futures
futures_list <- list()
for (c in 1:n_chains) {
  futures_list[[c]] <- future({
    run_block_adaptive_chain(Y, n_iter, inits_list[[c]], c, progress_dir)
  }, seed = TRUE)
}

# Monitor Progress
cat("   Monitoring chains (Acceptance should settle around 20-30%)...\n")
while (!all(resolved(futures_list))) {
  status_texts <- character(n_chains)
  for (c in 1:n_chains) {
    fpath <- file.path(progress_dir, paste0("status_", c, ".txt"))
    if (file.exists(fpath)) {
      info <- suppressWarnings(try(readLines(fpath, n = 1), silent = TRUE))
      if (inherits(info, "try-error") || length(info) == 0) {
        status_texts[c] <- "Init..."
      } else {
        parts <- strsplit(info, "\\|")[[1]]
        status_texts[c] <- sprintf("C%d: %3s%% (Acc: %4s%%)", c, parts[1], parts[2])
      }
    } else {
      status_texts[c] <- sprintf("C%d: Pending...", c)
    }
  }
  cat("\r", paste(status_texts, collapse = " | "))
  flush.console()
  Sys.sleep(0.5)
}
cat("\nDone.\n")

# -----------------------------------------------------------
# 5. Analysis & Diagnostics
# -----------------------------------------------------------
chains_list <- value(futures_list)
mcmc_obj <- mcmc.list(chains_list)
mcmc_clean <- window(mcmc_obj, start = burn_in + 1, thin = 10)

cat("\n--- DIAGNOSTICS ---\n")
gelman_res <- gelman.diag(mcmc_clean)
print(gelman_res)

eff_size <- effectiveSize(mcmc_clean)
cat("\nEffective Sample Sizes (Theta):", round(eff_size["theta"]), "\n")

# -----------------------------------------------------------
# 6. Plotting
# -----------------------------------------------------------
par(mfrow = c(2, 1))
plot(mcmc_clean[, "theta"], main = "Traceplot: Theta", auto.layout = FALSE)
par(mfrow = c(1, 1))

dens_est <- density(as.matrix(mcmc_clean)[, "theta"])
plot(dens_est,
  main = "Posterior Density of Theta (Full Likelihood + Adaptive)", 
  lwd = 2, col = "blue",
  xlab = expression(theta)
)
abline(v = theta_true, col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)

# Posterior Summary
cat("\nPosterior Mean Theta:", round(mean(as.matrix(mcmc_clean[, "theta"])), 3), "\n")
cat("Posterior Median Theta:", round(median(as.matrix(mcmc_clean[, "theta"])), 3), "\n")


plot(mcmc_clean)
