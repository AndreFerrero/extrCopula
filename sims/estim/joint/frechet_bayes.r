library(copula)
library(coda)
library(mvtnorm)
library(future)
library(future.apply)

set.seed(123)

# -----------------------------------------------------------
# 1. Simulate Data (Fréchet Margins)
# -----------------------------------------------------------
n <- 150
theta_true <- 3

# Fréchet Parameters
# shape (alpha), location (mu), scale (sigma)
alpha_true <- 3   # Shape (Tail index)
mu_true    <- 10  # Location (Lower bound)
sigma_true <- 5   # Scale

cat("1. Simulating Data (Fréchet Margins, n =", n, ")...\n")

# A. Dependence
cop <- gumbelCopula(param = theta_true, dim = n)
U <- as.numeric(rCopula(1, cop))

# B. Margins (Fréchet Quantile Function)
# X = mu + sigma * (-log(U))^(-1/alpha)
X <- mu_true + sigma_true * (-log(U))^(-1/alpha_true)

# -----------------------------------------------------------
# 2. Model Definitions (Fréchet + Gumbel)
# -----------------------------------------------------------

# Custom Density for Shifted Fréchet
d_frechet <- function(x, mu, sigma, alpha) {
  # Support check handled in posterior, but safety here:
  z <- (x - mu) / sigma
  if(any(z <= 0)) return(rep(0, length(x)))
  
  # PDF: (alpha/sigma) * z^(-1-alpha) * exp(-z^-alpha)
  (alpha / sigma) * (z^(-1 - alpha)) * exp(-(z^(-alpha)))
}

# Custom CDF for Shifted Fréchet
p_frechet <- function(x, mu, sigma, alpha) {
  z <- (x - mu) / sigma
  # If x <= mu, probability is 0
  p <- exp(-(z^(-alpha)))
  p[x <= mu] <- 0
  return(p)
}

log_posterior <- function(param_vec, data) {
  # Unpack parameters
  theta <- param_vec[1]
  mu    <- param_vec[2]
  sigma <- param_vec[3]
  alpha <- param_vec[4]

  # --- A. CONSTRAINTS (CRITICAL) ---
  
  # 1. Copula constraints
  if (theta <= 1.01 || theta > 20) return(-Inf)
  
  # 2. Fréchet Constraints
  if (sigma <= 0.01) return(-Inf)
  if (alpha <= 0.1 || alpha > 20) return(-Inf)
  
  # 3. SUPPORT CONSTRAINT: 
  # Fréchet is only defined for x > mu. 
  # If any data point is <= mu, likelihood is 0 (log is -Inf).
  if (mu >= min(data)) return(-Inf)

  # --- B. PRIORS ---
  
  lp_theta <- dgamma(theta - 1, shape = 2, rate = 1, log = TRUE)
  
  # Location (mu): Normal, but effectively truncated by the data min
  lp_mu <- dnorm(mu, mean = mean(data) - sd(data), sd = 10, log = TRUE)
  
  # Scale (sigma): Boundary Avoiding LogNormal
  lp_sigma <- dlnorm(sigma, meanlog = log(sd(data)/2), sdlog = 0.5, log = TRUE)
  
  # Shape (alpha): Gamma (favoring values around 2-5 typical for extremes)
  lp_alpha <- dgamma(alpha, shape = 3, rate = 1, log = TRUE)

  # --- C. MARGINAL LIKELIHOOD ---
  dens_vals <- d_frechet(data, mu, sigma, alpha)
  
  # Safety for log(0)
  if (any(dens_vals <= 0)) return(-Inf)
  ll_margin <- sum(log(dens_vals))

  # --- D. COPULA LIKELIHOOD ---
  u_hat <- p_frechet(data, mu, sigma, alpha)
  u_hat <- pmin(pmax(u_hat, 1e-8), 1 - 1e-8)

  ld_cop <- dCopula(u_hat, 
                    copula = gumbelCopula(theta, dim = length(data)), 
                    log = TRUE)
  
  if (is.na(ld_cop) || !is.finite(ld_cop)) return(-Inf)

  return(lp_theta + lp_mu + lp_sigma + lp_alpha + ll_margin + ld_cop)
}

# -----------------------------------------------------------
# 3. Worker: Block Adaptive Metropolis
# -----------------------------------------------------------
run_worker <- function(data, n_iter, init_vals, chain_id, progress_dir) {
  library(copula); library(mvtnorm); library(coda)
  
  p_names <- c("theta", "mu", "sigma", "alpha")
  n_par <- length(p_names)
  chain <- matrix(NA, nrow = n_iter, ncol = n_par)
  colnames(chain) <- p_names

  curr_par <- init_vals
  
  # Robust Start
  curr_lp <- log_posterior(curr_par, data)
  if(!is.finite(curr_lp)) {
     stop("Worker failed: Initial values yield -Inf posterior (Likely mu >= min(x)).")
  }

  # Adaptive Settings
  cov_mat <- diag(rep(0.01, n_par)) 
  # Optimal Scaling for d=4
  sd_scale <- (2.38^2) / n_par 
  eps <- 1e-6 * diag(n_par)
  
  adapt_start <- 500
  total_accepts <- 0
  status_file <- file.path(progress_dir, paste0("status_", chain_id, ".txt"))

  for (i in 1:n_iter) {
    # Block Proposal
    prop_val <- rmvnorm(1, mean = curr_par, sigma = cov_mat)[1, ]
    
    prop_lp <- log_posterior(prop_val, data)
    ratio <- prop_lp - curr_lp

    if (is.finite(ratio) && log(runif(1)) < ratio) {
      curr_par <- prop_val
      curr_lp <- prop_lp
      total_accepts <- total_accepts + 1
    }
    chain[i, ] <- curr_par

    # Adaptation
    if (i > adapt_start && i %% 50 == 0) {
      emp_cov <- cov(chain[1:i, ])
      cov_mat <- (sd_scale * emp_cov) + eps
    }

    # Dashboard
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

n_iter <- 12000
burn_in <- 3000

inits_list <- list()
for (c in 1:n_chains) {
  # Careful Initialization for mu: Must be < min(X)
  data_min <- min(X)
  
  inits_list[[c]] <- c(
    theta = runif(1, 2, 4), 
    mu    = data_min - runif(1, 0.5, 2), # Ensure below min
    sigma = runif(1, 2, 6),
    alpha = runif(1, 2, 4)
  )
}

futures_list <- list()
for (c in 1:n_chains) {
  futures_list[[c]] <- future({
    run_worker(X, n_iter, inits_list[[c]], c, progress_dir)
  }, seed = TRUE)
}

cat("MCMC Started (Fréchet Margins)...\n")

# Robust Monitor
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
res_dir <- here("sims", "estim", "joint", "res")

save(futures_list, file = here(res_dir, "frechet_bayes_chains.Rdata"))
load(here(res_dir, "frechet_bayes_chains.Rdata"))

chains_list <- value(futures_list)
mcmc_obj <- mcmc.list(chains_list)
mcmc_clean <- window(mcmc_obj, start = burn_in + 1, thin = 5)

cat("\n--- DIAGNOSTICS ---\n")
print(gelman.diag(mcmc_clean))
print(effectiveSize(mcmc_clean))

# Visuals
par(mfrow = c(2, 2))
plot(mcmc_clean[, "theta"], main="Theta", density=FALSE, auto.layout=FALSE)
plot(mcmc_clean[, "mu"], main="Mu (Loc)", density=FALSE, auto.layout=FALSE)
plot(mcmc_clean[, "sigma"], main="Sigma (Scale)", density=FALSE, auto.layout=FALSE)
plot(mcmc_clean[, "alpha"], main="Alpha (Shape)", density=FALSE, auto.layout=FALSE)

par(mfrow=c(1,1))
dens_est <- density(as.matrix(mcmc_clean)[, "theta"])
plot(dens_est, main = "Posterior Theta (Fréchet Margin)", lwd = 2, col = "blue", xlab=expression(theta))
abline(v = theta_true, col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)

cat("\nTrue Theta:", theta_true)
cat("\nPosterior Mean:", round(mean(as.matrix(mcmc_clean[,"theta"])), 3))


# Blocking analysis: should I block some parameters
library(coda)
library(MASS)   # kde2d
library(psych)  # pairs.panels (optional, nice viz)

# convert to matrix of samples (chains combined)
mat_all <- as.matrix(mcmc_clean)        # or as.matrix(mcmc_clean) if you want post-burnin
pnames <- colnames(mat_all)

# 1.A: Posterior correlation matrix
cor_mat <- cor(mat_all)
print(round(cor_mat, 3))

# 1.B: Spearman (rank) correlations (robust to nonlinearity)
spearman_mat <- cor(mat_all, method = "spearman")
print(round(spearman_mat, 3))

# 1.C: Pairs plot with densities and smoothing (nice quick view)
pairs.panels(mat_all[, c("sigma","alpha","mu","theta")], ellipses = FALSE, lm = TRUE)

# 1.D: Bivariate density contour for sigma vs alpha
z <- MASS::kde2d(mat_all[,"sigma"], mat_all[,"alpha"], n = 80)
plot(mat_all[,"sigma"], mat_all[,"alpha"], pch = 20, cex = 0.4,
     xlab="sigma", ylab="alpha", main="sigma vs alpha (posterior)")
contour(z, add = TRUE)

z <- MASS::kde2d(mat_all[,"theta"], mat_all[,"alpha"], n = 80)
plot(mat_all[,"theta"], mat_all[,"alpha"], pch = 20, cex = 0.4,
     xlab="theta", ylab="alpha", main="theta vs alpha (posterior)")
contour(z, add = TRUE)
