library(copula)
library(coda)
library(mvtnorm)
library(future)
library(future.apply)

set.seed(123)

# -----------------------------------------------------------
# 1. Simulate High-Dimensional Data (n=300)
# -----------------------------------------------------------
n_dim <- 300
theta_true <- 4
mu_true <- 0
sigma_true <- 1

cat("Simulating data (dim =", n_dim, ")...\n")
cop <- gumbelCopula(param = theta_true, dim = n_dim)
U <- as.numeric(rCopula(1, cop))
X <- qlnorm(U, meanlog = mu_true, sdlog = sigma_true)
Y <- log(X)

# -----------------------------------------------------------
# 2. Pairwise Likelihood Helper
# -----------------------------------------------------------
# Pre-calculate indices.
# Full pairwise for n=300 is ~45,000 pairs. This is okay, but
# can be slow inside MCMC. We subsample for speed here.

all_pairs <- t(combn(n_dim, 2))
n_pairs_total <- nrow(all_pairs)

# Subsampling strategy: Keep max 2000 pairs for speed
n_keep <- 3000

if (n_pairs_total > n_keep) {
  set.seed(42)
  keep_idx <- sample(n_pairs_total, n_keep)
  pair_indices <- all_pairs[keep_idx, ]
  cat("Subsampling:", n_keep, "pairs out of", n_pairs_total, "\n")
} else {
  pair_indices <- all_pairs
}

calc_pairwise_ll <- function(u_vec, theta) {
  # 1. Extract the pairs from the vector u_vec
  # resulting matrix has 2 columns
  u_pairs <- cbind(u_vec[pair_indices[, 1]], u_vec[pair_indices[, 2]])

  # 2. Vectorized calculation of bivariate densities
  # dCopula for dim=2 is fast/optimized in C
  d_vals <- dCopula(u_pairs, copula = gumbelCopula(theta, dim = 2))

  # 3. Handle numerical issues
  if (any(is.na(d_vals)) || any(d_vals <= 0)) {
    return(-1e10)
  }

  # 4. Sum log-densities
  return(sum(log(d_vals)))
}

# -----------------------------------------------------------
# 3. Bayesian Model (Mixture Margin + Pairwise Copula)
# -----------------------------------------------------------
d_mix <- function(y, w, m1, s1, m2, s2) {
  w * dnorm(y, m1, s1) + (1 - w) * dnorm(y, m2, s2)
}

p_mix <- function(y, w, m1, s1, m2, s2) {
  w * pnorm(y, m1, s1) + (1 - w) * pnorm(y, m2, s2)
}

log_posterior <- function(param_vec, data_y) {
  theta <- param_vec[1]
  w <- param_vec[2]
  mu1 <- param_vec[3]
  s1 <- param_vec[4]
  mu2 <- param_vec[5]
  s2 <- param_vec[6]

  # --- PRIORS ---
  if (theta <= 1.01 || theta > 20) {
    return(-Inf)
  }
  lp_theta <- dgamma(theta - 1, shape = 2, rate = 1, log = TRUE)

  if (w <= 0.01 || w >= 0.99) {
    return(-Inf)
  }
  if (mu1 > mu2) {
    return(-Inf)
  }
  if (s1 <= 0.05 || s2 <= 0.05) {
    return(-Inf)
  }

  lp_s <- dgamma(s1, 3, 2, log = TRUE) + dgamma(s2, 3, 2, log = TRUE)
  lp_mu <- dnorm(mu1, 0, 5, log = TRUE) + dnorm(mu2, 0, 5, log = TRUE)

  # --- MARGINAL LIKELIHOOD ---
  dens_vals <- d_mix(data_y, w, mu1, s1, mu2, s2)
  if (any(dens_vals <= 1e-300)) {
    return(-Inf)
  }
  ll_margin <- sum(log(dens_vals))

  # --- TRANSFORM TO U ---
  u_hat <- p_mix(data_y, w, mu1, s1, mu2, s2)
  u_hat <- pmin(pmax(u_hat, 1e-6), 1 - 1e-6)

  # --- PAIRWISE COMPOSITE LIKELIHOOD ---
  ll_cop <- calc_pairwise_ll(u_hat, theta)

  return(lp_theta + lp_s + lp_mu + ll_margin + ll_cop)
}

# -----------------------------------------------------------
# 4. Parallel Worker
# -----------------------------------------------------------
run_pairwise_worker <- function(data, n_iter, init_vals, chain_id, progress_dir) {
  library(copula)
  library(mvtnorm)
  library(coda)

  p_names <- c("theta", "w", "mu1", "s1", "mu2", "s2")
  n_par <- length(p_names)
  chain <- matrix(NA, nrow = n_iter, ncol = n_par)

  colnames(chain) <- p_names

  curr_par <- init_vals
  curr_lp <- log_posterior(curr_par, data)

  prop_sd <- c(0.2, 0.05, 0.1, 0.1, 0.1, 0.1)
  batch_size <- 50
  accept_batch <- numeric(n_par)
  total_accepts <- 0
  status_file <- file.path(progress_dir, paste0("status_", chain_id, ".txt"))

  for (i in 1:n_iter) {
    for (j in 1:n_par) {
      prop_val <- curr_par
      prop_val[j] <- rnorm(1, curr_par[j], prop_sd[j])
      prop_lp <- log_posterior(prop_val, data)
      ratio <- prop_lp - curr_lp

      if (is.finite(ratio) && log(runif(1)) < ratio) {
        curr_par <- prop_val
        curr_lp <- prop_lp
        accept_batch[j] <- accept_batch[j] + 1
        total_accepts <- total_accepts + 1
      }
    }
    chain[i, ] <- curr_par

    # Adaptation
    if (i <= n_iter / 2 && i %% batch_size == 0) {
      acc_rates <- accept_batch / batch_size
      for (k in 1:n_par) {
        if (acc_rates[k] < 0.2) {
          prop_sd[k] <- prop_sd[k] * 0.8
        } else if (acc_rates[k] > 0.3) prop_sd[k] <- prop_sd[k] * 1.2
      }
      accept_batch <- numeric(n_par)
    }

    # Dashboard
    if (i %% 50 == 0 || i == n_iter) {
      rate <- (total_accepts / (i * n_par)) * 100
      pct <- as.integer((i / n_iter) * 100)
      tryCatch(
        {
          writeLines(sprintf("%d|%.1f", pct, rate), status_file)
        },
        error = function(e) NULL
      )
    }
  }
  return(mcmc(chain))
}

# -----------------------------------------------------------
# 5. Execution
# -----------------------------------------------------------
plan(multisession, workers = 3)
progress_dir <- tempdir()
file.remove(list.files(progress_dir, pattern = "status_", full.names = TRUE))

cat("Launching Pairwise MCMC...\n")
n_chains <- 3
n_iter <- 5000
burn_in <- 1000

inits_list <- list()
for (c in 1:n_chains) {
  inits_list[[c]] <- c(
    theta = runif(1, 2, 5), w = 0.5,
    mu1 = mean(Y) - 0.2, s1 = sd(Y), mu2 = mean(Y) + 0.2, s2 = sd(Y)
  )
}

futures_list <- list()
for (c in 1:n_chains) {
  futures_list[[c]] <- future(
    {
      run_pairwise_worker(Y, n_iter, inits_list[[c]], c, progress_dir)
    },
    seed = TRUE
  )
}

# Monitor
while (!all(resolved(futures_list))) {
  status_texts <- character(n_chains)
  for (c in 1:n_chains) {
    fpath <- file.path(progress_dir, paste0("status_", c, ".txt"))
    if (file.exists(fpath)) {
      info <- suppressWarnings(try(readLines(fpath, n = 1), silent = TRUE))
      if (!inherits(info, "try-error") && length(info) > 0) {
        parts <- strsplit(info, "\\|")[[1]]
        status_texts[c] <- sprintf("C%d: %3s%% (Acc: %4s%%)", c, parts[1], parts[2])
      } else {
        status_texts[c] <- "Init..."
      }
    } else {
      status_texts[c] <- "Pending..."
    }
  }
  cat("\r", paste(status_texts, collapse = " | "))
  flush.console()
  Sys.sleep(0.5)
}
cat("\nDone.\n")

# -----------------------------------------------------------
# 6. Results
# -----------------------------------------------------------
chains_list <- value(futures_list)
mcmc_clean <- window(mcmc.list(chains_list), start = burn_in + 1)
summary_stats <- summary(mcmc_clean)$quantiles
cat("\n--- PAIRWISE RESULTS ---\n")
print(round(summary_stats[, c(1, 3, 5)], 3))

dens_theta <- density(as.matrix(mcmc_clean)[, 1])
plot(dens_theta, main = "Pairwise Posterior", col = "blue", lwd = 2)
abline(v = theta_true, col = "red", lwd = 2, lty = 2)
