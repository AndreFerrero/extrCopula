library(copula)
library(coda)
library(mvtnorm)
library(future)
library(future.apply)

set.seed(123)

# -----------------------------------------------------------
# 1. Simulate Data (Direct Student-t)
# -----------------------------------------------------------
n <- 150
theta_true <- 3

# Student-t Parameters
mu_true <- 10
sigma_true <- 2
nu_true <- 5

cat("1. Simulating Data (n =", n, ")...\n")

# A. Dependence
cop <- gumbelCopula(param = theta_true, dim = n)
U <- as.numeric(rCopula(1, cop))

# B. Margins (Transform U to X directly via Quantile Function)
X <- mu_true + sigma_true * qt(U, df = nu_true)


# -----------------------------------------------------------
# 2. Model Definitions
# -----------------------------------------------------------

# Helper for Location-Scale t-density
d_t_ls <- function(x, m, s, df) {
  dt((x - m) / s, df = df) / s
}

# Helper for Location-Scale t-CDF
p_t_ls <- function(x, m, s, df) {
  pt((x - m) / s, df = df)
}

log_posterior <- function(param_vec, data) {
  # Unpack parameters
  theta <- param_vec[1]
  mu <- param_vec[2]
  sigma <- param_vec[3]
  nu <- param_vec[4]

  # --- A. PRIORS ---

  # Copula: Gamma(2,1) shifted
  if (theta <= 1.01 || theta > 20) {
    return(-Inf)
  }
  lp_theta <- dgamma(theta - 1, shape = 2, rate = 1, log = TRUE)

  # Location (mu): Vague Normal
  lp_mu <- dnorm(mu, mean = mean(data), sd = 10, log = TRUE)

  # Scale (sigma): Boundary Avoiding (Log-Normal)
  # Prevents sigma -> 0 singularity
  if (sigma <= 0.01) {
    return(-Inf)
  }
  lp_sigma <- dlnorm(sigma, meanlog = log(sd(data)), sdlog = 0.5, log = TRUE)

  # Degrees of Freedom (nu): Penalized
  # Restrict to > 2.1 for finite variance, upper bound for stability
  if (nu <= 2.1 || nu > 100) {
    return(-Inf)
  }
  lp_nu <- dexp(nu - 2, rate = 0.1, log = TRUE)

  # --- B. MARGINAL LIKELIHOOD ---
  # Sum of log-densities
  dens_vals <- d_t_ls(data, mu, sigma, nu)

  # Safety check
  if (any(dens_vals <= 0)) {
    return(-Inf)
  }
  ll_margin <- sum(log(dens_vals))

  # --- C. COPULA LIKELIHOOD (Full Density) ---
  # 1. Transform Data to U
  u_hat <- p_t_ls(data, mu, sigma, nu)

  # 2. Clamp for numerical safety (avoid exactly 0 or 1)
  u_hat <- pmin(pmax(u_hat, 1e-8), 1 - 1e-8)

  # 3. Calculate Full Likelihood
  ld_cop <- dCopula(u_hat,
    copula = gumbelCopula(theta, dim = length(data)),
    log = TRUE
  )

  # 4. Check for numerical failure (underflow/overflow)
  if (is.na(ld_cop) || !is.finite(ld_cop)) {
    return(-Inf)
  }

  # Total Posterior
  return(lp_theta + lp_mu + lp_sigma + lp_nu + ll_margin + ld_cop)
}

# -----------------------------------------------------------
# 3. Worker: Block Adaptive Metropolis
# -----------------------------------------------------------
run_worker <- function(data, n_iter, init_vals, chain_id, progress_dir) {
  library(copula)
  library(mvtnorm)
  library(coda)

  p_names <- c("theta", "mu", "sigma", "nu")
  n_par <- length(p_names)
  chain <- matrix(NA, nrow = n_iter, ncol = n_par)
  colnames(chain) <- p_names

  curr_par <- init_vals

  # Robust Start: Ensure we don't start at -Inf
  curr_lp <- log_posterior(curr_par, data)
  if (!is.finite(curr_lp)) {
    stop("Worker failed: Initial values yield -Inf posterior. Try better inits.")
  }

  # --- ADAPTIVE METROPOLIS SETTINGS ---
  # Initial small covariance
  cov_mat <- diag(rep(0.01, n_par))

  # Optimal Scaling for d=4 (Haario et al)
  sd_scale <- (2.38^2) / n_par
  eps <- 1e-6 * diag(n_par) # Regularization

  adapt_start <- 500
  total_accepts <- 0
  status_file <- file.path(progress_dir, paste0("status_", chain_id, ".txt"))

  for (i in 1:n_iter) {
    # 1. BLOCK PROPOSAL
    prop_val <- rmvnorm(1, mean = curr_par, sigma = cov_mat)[1, ]

    prop_lp <- log_posterior(prop_val, data)
    ratio <- prop_lp - curr_lp

    # 2. ACCEPT/REJECT
    if (is.finite(ratio) && log(runif(1)) < ratio) {
      curr_par <- prop_val
      curr_lp <- prop_lp
      total_accepts <- total_accepts + 1
    }
    chain[i, ] <- curr_par

    # 3. ADAPTATION
    if (i > adapt_start && i %% 50 == 0) {
      emp_cov <- cov(chain[1:i, ])
      cov_mat <- (sd_scale * emp_cov) + eps
    }

    # 4. DASHBOARD
    if (i %% 50 == 0 || i == n_iter) {
      rate <- (total_accepts / i) * 100
      pct <- as.integer((i / n_iter) * 100)
      try(
        {
          writeLines(sprintf("%d|%.1f", pct, rate), status_file)
        },
        silent = TRUE
      )
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

inits_list <- list()
for (c in 1:n_chains) {
  inits_list[[c]] <- c(
    theta = runif(1, 2, 4),
    mu = mean(X),
    sigma = sd(X),
    nu = 10
  )
}

futures_list <- list()
for (c in 1:n_chains) {
  futures_list[[c]] <- future(
    {
      run_worker(X, n_iter, inits_list[[c]], c, progress_dir)
    },
    seed = TRUE
  )
}

cat("MCMC Started (Full Likelihood + Student-t Margin)...\n")

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
      status_texts[c] <- "Init..."
    }
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
load(here(res_dir, "t_bayes_chains.Rdata"))

chains_list <- value(futures_list)
mcmc_obj <- mcmc.list(chains_list)

burn_in <- n_iter/2
mcmc_clean <- window(mcmc_obj, start = burn_in + 1, thin = 5)

cat("\n--- DIAGNOSTICS ---\n")
print(gelman.diag(mcmc_clean))
print(effectiveSize(mcmc_clean))

# Visuals
par(mfrow = c(2, 2))
plot(mcmc_clean[, "theta"], main = "Theta", density = FALSE, auto.layout = FALSE)
plot(mcmc_clean[, "sigma"], main = "Sigma", density = FALSE, auto.layout = FALSE)
plot(mcmc_clean[, "mu"], main = "Mu", density = FALSE, auto.layout = FALSE)
plot(mcmc_clean[, "nu"], main = "Nu", density = FALSE, auto.layout = FALSE)
par(mfrow = c(1, 1))

paste("True Theta:", theta_true)
paste("Posterior Mean:", round(mean(as.matrix(mcmc_clean[, "theta"])), 3))

theta_dens_est <- density(as.matrix(mcmc_clean)[, "theta"])
plot(theta_dens_est,
  main = "Posterior Density of Theta", lwd = 2, col = "blue",
  xlab = expression(theta)
)
abline(v = theta_true, col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)

nu_dens_est <- density(as.matrix(mcmc_clean)[, "nu"])
plot(nu_dens_est,
  main = "Posterior Density of Nu", lwd = 2, col = "blue",
  xlab = expression(nu)
)
abline(v = nu_true, col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)

sigma_dens_est <- density(as.matrix(mcmc_clean)[, "sigma"])
plot(sigma_dens_est,
  main = "Posterior Density of Sigma", lwd = 2, col = "blue",
  xlab = expression(sigma)
)
abline(v = sigma_true, col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)