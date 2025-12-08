library(copula)
library(coda)
library(mvtnorm)
library(future)
library(future.apply)
library(sn) # CRITICAL: For Skew-Normal distribution

set.seed(123)

# -----------------------------------------------------------
# 1. Simulate Data (Same as before)
# -----------------------------------------------------------
n <- 150
theta_true <- 3
mu_true <- 0
sigma_true <- 1

cop <- gumbelCopula(param = theta_true, dim = n)
U <- as.numeric(rCopula(1, cop))
X <- qlnorm(U, meanlog = mu_true, sdlog = sigma_true)
Y <- log(X)

# -----------------------------------------------------------
# 2. Skew-Normal Model Definitions
# -----------------------------------------------------------
# We replace the 5-parameter Mixture with a 3-parameter Skew-Normal
# Params: xi (loc), omega (scale), alpha (shape/skew)

log_posterior <- function(param_vec, data_y, scaling_k) {
  # Unpack
  theta <- param_vec[1]
  xi    <- param_vec[2]
  omega <- param_vec[3]
  alpha <- param_vec[4]

  # --- PRIORS ---
  if (theta <= 1.01 || theta > 20) return(-Inf)
  lp_theta <- dgamma(theta - 1, shape = 2, rate = 1, log = TRUE)

  # Prior on Omega (Scale): Strict non-zero prior (Inverse Gamma-like)
  # This prevents singularity just like in the mixture model
  if (omega <= 0.05) return(-Inf)
  lp_omega <- dgamma(omega, shape = 2, rate = 2, log = TRUE) 
  
  # Prior on Xi (Location): Normal
  lp_xi <- dnorm(xi, mean=0, sd=5, log=TRUE)
  
  # Prior on Alpha (Skewness): Normal (centered at 0 = Normal distribution)
  lp_alpha <- dnorm(alpha, mean=0, sd=5, log=TRUE)

  # --- MARGINAL LIKELIHOOD ---
  # dsn handles the skew-normal density
  # log=TRUE is more stable
  dens_vals <- dsn(data_y, xi=xi, omega=omega, alpha=alpha, log=TRUE)
  ll_margin <- sum(dens_vals)

  # --- COPULA LIKELIHOOD ---
  # Transform data to Uniforms using psn (CDF)
  u_hat <- psn(data_y, xi=xi, omega=omega, alpha=alpha)
  
  # Clamp for safety
  u_hat <- pmin(pmax(u_hat, 1e-7), 1 - 1e-7)

  # Calculate Pairwise Composite Likelihood
  # (Using the pairwise method is safer than full likelihood for stability)
  # Recalculating pairwise indices locally for this function scope
  # (In production, pass these in or ensure global env visibility)
  
  # ** NOTE: For this example, I will use the Full Likelihood since n=150 **
  # ** If n > 200, switch back to calc_pairwise_ll **
  
  ld_cop <- dCopula(u_hat, 
                    copula = gumbelCopula(theta, dim = length(data_y)), 
                    log = TRUE)
  
  if (is.na(ld_cop) || !is.finite(ld_cop)) return(-Inf)

  # Note: No 'scaling_k' needed for Full Likelihood. 
  # If using Composite Pairwise, multiply ld_cop by k.
  
  return(lp_theta + lp_omega + lp_xi + lp_alpha + ll_margin + ld_cop)
}

# -----------------------------------------------------------
# 3. Parallel Worker
# -----------------------------------------------------------
run_worker_sn <- function(data, n_iter, init_vals, chain_id, progress_dir) {
  library(copula); library(mvtnorm); library(coda); library(sn)
  
  p_names <- c("theta", "xi", "omega", "alpha")
  n_par <- length(p_names)
  chain <- matrix(NA, nrow = n_iter, ncol = n_par)
  colnames(chain) <- p_names

  curr_par <- init_vals
  
  # Initialization Check
  curr_lp <- log_posterior(curr_par, data, 1.0)
  if(!is.finite(curr_lp)) stop("Invalid start")

  # --- ADAPTIVE METROPOLIS ---
  # Initial covariance
  cov_mat <- diag(rep(0.01, n_par)) 
  # Reduced step size slightly to fix your "sticky" traceplot
  # (1.8^2) / d is safer than (2.38^2) / d for non-gaussian targets
  sd_scale <- (1.8^2) / n_par 
  eps <- 1e-6 * diag(n_par)
  
  adapt_start <- 500
  total_accepts <- 0
  status_file <- file.path(progress_dir, paste0("status_", chain_id, ".txt"))

  for (i in 1:n_iter) {
    prop_val <- rmvnorm(1, mean = curr_par, sigma = cov_mat)[1, ]
    prop_lp <- log_posterior(prop_val, data, 1.0)
    ratio <- prop_lp - curr_lp

    if (is.finite(ratio) && log(runif(1)) < ratio) {
      curr_par <- prop_val
      curr_lp <- prop_lp
      total_accepts <- total_accepts + 1
    }
    chain[i, ] <- curr_par

    if (i > adapt_start && i %% 50 == 0) {
      emp_cov <- cov(chain[1:i, ])
      cov_mat <- (sd_scale * emp_cov) + eps
    }

    if (i %% 50 == 0 || i == n_iter) {
      rate <- (total_accepts / i) * 100
      pct <- as.integer((i / n_iter) * 100)
      try({ writeLines(sprintf("%d|%.1f", pct, rate), status_file) }, silent=TRUE)
    }
  }
  return(mcmc(chain))
}

# -----------------------------------------------------------
# 4. Execution
# -----------------------------------------------------------
n_chains <- 3
plan(multisession, workers = n_chains)
progress_dir <- tempdir()
file.remove(list.files(progress_dir, pattern = "status_", full.names = TRUE))

# We can run FEWER iterations because mixing will be much faster
n_iter <- 8000 
burn_in <- 1000

inits_list <- list()
for (c in 1:n_chains) {
  inits_list[[c]] <- c(
    theta = runif(1, 2, 4), 
    xi = mean(Y), 
    omega = sd(Y), 
    alpha = 0 # Start with assumption of Normality
  )
}

futures_list <- list()
for (c in 1:n_chains) {
  futures_list[[c]] <- future({
    run_worker_sn(Y, n_iter, inits_list[[c]], c, progress_dir)
  }, seed = TRUE)
}

cat("MCMC Started (Skew-Normal Model)...\n")
while (!all(resolved(futures_list))) {
  # ... (Same dashboard code as before) ...
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
    } else { status_texts[c] <- "Init..." }
  }
  cat("\r", paste(status_texts, collapse = " | "))
  flush.console()
  Sys.sleep(0.5)
}
cat("\nDone.\n")

# -----------------------------------------------------------
# 5. Analysis
# -----------------------------------------------------------
chains_list <- value(futures_list)
mcmc_obj <- mcmc.list(chains_list)
mcmc_clean <- window(mcmc_obj, start = burn_in + 1, thin = 5)

cat("\n--- DIAGNOSTICS (Skew-Normal) ---\n")
print(gelman.diag(mcmc_clean))
print(effectiveSize(mcmc_clean))

par(mfrow = c(2, 2))
plot(mcmc_clean) # Check mixing for ALL parameters

dens_est <- density(as.matrix(mcmc_clean)[, "theta"])
plot(dens_est,
  main = "Posterior Density of Theta", lwd = 2, col = "blue",
  xlab = expression(theta)
)
abline(v = theta_true, col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)