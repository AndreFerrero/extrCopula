library(copula)
library(evd)
library(coda)
library(mvtnorm)
library(future)
library(future.apply)
library(MASS) # for kde2d
library(progressr)

set.seed(123)

# -----------------------------------------------------------
# 1. Simulate Data (Fréchet Margins with Mu = 0)
# -----------------------------------------------------------
n <- 200
theta_true <- 1.5

# Fréchet Parameters
alpha_true <- 3   # Shape (Tail index)
mu_true    <- 0   # FIXED TO 0 (Standard Fréchet lower bound)
sigma_true <- 5   # Scale

cat("1. Simulating Data (Fréchet Margins, mu=0, n =", n, ")...\n")

# A. Dependence
cop <- gumbelCopula(param = theta_true, dim = n)
U <- as.numeric(rCopula(1, cop))

# B. Margins
# X = mu + sigma * (-log(U))^(-1/alpha)
X <- mu_true + sigma_true * (-log(U))^(-1 / alpha_true)

# -----------------------------------------------------------
# 2. Transformed Log-Posterior (Strictly Necessary Change)
# -----------------------------------------------------------
# We sample in transformed space to avoid boundary degeneracy.
# params_trans: [log(theta-1), log(sigma), log(alpha)]

log_posterior_transformed <- function(param_trans, data) {
  
  # 1. BACK-TRANSFORM to Natural Scale
  theta <- exp(param_trans[1]) + 1
  sigma <- exp(param_trans[2])
  alpha <- exp(param_trans[3])
  
  # Fixed location
  mu_fixed <- mu_true

  # 2. JACOBIAN ADJUSTMENT
  # log|J| = sum(log(derivs)) = sum(log(val)) for log-transform
  # derivative of (exp(y)+1) is exp(y) = theta-1
  # derivative of exp(y) is exp(y) = sigma (or alpha)
  # Simply: sum of the transformed params (since p_trans[1] = log(theta-1))
  log_jacobian <- param_trans[1] + param_trans[2] + param_trans[3]

  # 3. PRIORS (Defined on Natural Scale)
  lp_theta <- dgamma(theta - 1, shape = 2, rate = 1, log = TRUE)
  lp_sigma <- dlnorm(sigma, meanlog = log(sd(data)/2), sdlog = 1.0, log = TRUE)
  lp_alpha <- dgamma(alpha, shape = 3, rate = 1, log = TRUE)

  # 4. MARGINAL LIKELIHOOD
  if (sigma < 1e-10 || alpha < 1e-10) return(-Inf)
  
  dens_vals <- dfrechet(data, loc = mu_fixed, scale = sigma, shape = alpha)
  if (any(dens_vals <= 0)) return(-Inf)
  ll_margin <- sum(log(dens_vals))

  # 5. COPULA LIKELIHOOD
  u_hat <- pfrechet(data, loc = mu_fixed, scale = sigma, shape = alpha)
  u_hat <- pmin(pmax(u_hat, 1e-9), 1 - 1e-9)

  ld_cop <- dCopula(u_hat, copula = gumbelCopula(theta, dim = length(data)), log = TRUE)
  
  if (!is.finite(ld_cop)) return(-Inf)

  return(lp_theta + lp_sigma + lp_alpha + ll_margin + ld_cop + log_jacobian)
}

# -----------------------------------------------------------
# 3. Worker: Adaptive Metropolis on Transformed Space
# -----------------------------------------------------------
run_worker_trans <- function(data, n_iter, init_trans, chain_id, progress_dir) {
  library(copula)
  library(mvtnorm)
  library(coda)
  library(evd)

  p_names <- c("log_theta_m1", "log_sigma", "log_alpha")
  n_par <- length(p_names)
  
  # Matrix to store Transformed samples
  chain_trans <- matrix(NA, nrow = n_iter, ncol = n_par)

  curr_trans <- init_trans
  
  # Check init
  curr_lp <- log_posterior_transformed(curr_trans, data)
  if (!is.finite(curr_lp)) {
    stop("Worker failed: Initial values yield -Inf posterior.")
  }

  # Adaptive Settings
  cov_mat <- diag(rep(0.1, n_par)) # Small init variance
  sd_scale <- (2.38^2) / n_par
  eps <- 1e-6 * diag(n_par)

  adapt_start <- 500
  total_accepts <- 0
  status_file <- file.path(progress_dir, paste0("status_", chain_id, ".txt"))

  for (i in 1:n_iter) {
    # Proposal in Transformed Space (Unconstrained)
    prop_trans <- rmvnorm(1, mean = curr_trans, sigma = cov_mat)[1, ]

    prop_lp <- log_posterior_transformed(prop_trans, data)
    ratio <- prop_lp - curr_lp

    if (is.finite(ratio) && log(runif(1)) < ratio) {
      curr_trans <- prop_trans
      curr_lp <- prop_lp
      total_accepts <- total_accepts + 1
    }
    chain_trans[i, ] <- curr_trans

    # Adaptation
    if (i > adapt_start && i %% 50 == 0) {
      emp_cov <- cov(chain_trans[1:i, ])
      cov_mat <- (sd_scale * emp_cov) + eps
    }

    # Reporting
    if (i %% 500 == 0) {
      try({
        writeLines(sprintf("%d|%.1f", as.integer((i/n_iter)*100), (total_accepts/i)*100), status_file)
      }, silent = TRUE)
    }
  }
  return(chain_trans) # Return transformed chain
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
  # Init values on Natural Scale
  th_init <- runif(1, 2, 4)
  si_init <- runif(1, 3, 7)
  al_init <- runif(1, 2, 4)

  # Transform Inits
  inits_list[[c]] <- c(
    log(th_init - 1),
    log(si_init),
    log(al_init)
  )
}

futures_list <- list()
for (c in 1:n_chains) {
    futures_list[[c]] <- future({
      res <- run_worker_trans(X, n_iter, inits_list[[c]], c, progress_dir)
      cat(sprintf("Chain %d done\n", c))
      res
    }, seed = TRUE
  )
}

paste("MCMC Started (Correct Log-Normal Specification)...\n")

# Robust Monitoring Loop
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

chains_raw <- value(futures_list)

res_dir <- here("sims", "estim", "joint", "res")

# save(chains_raw, file = here(res_dir, "transf_frechet_bayes_chains.Rdata"))
load(here(res_dir, "transf_frechet_bayes_chains.Rdata"))

# -----------------------------------------------------------
# 5. Post-Processing: Back-Transform Chains
# -----------------------------------------------------------
chains_natural <- list()
for(c in 1:n_chains) {
  mat_trans <- chains_raw[[c]]
  mat_nat <- matrix(NA, nrow = nrow(mat_trans), ncol = 3)
  colnames(mat_nat) <- c("theta", "sigma", "alpha")
  
  mat_nat[, "theta"] <- exp(mat_trans[,1]) + 1
  mat_nat[, "sigma"] <- exp(mat_trans[,2])
  mat_nat[, "alpha"] <- exp(mat_trans[,3])
  
  chains_natural[[c]] <- mcmc(mat_nat)
}

mcmc_obj <- mcmc.list(chains_natural)
burn_in <- round(n_iter/3, 0)
mcmc_clean <- window(mcmc_obj, start = burn_in + 1, thin = 5)

# Analysis & Diagnostics
print(summary(mcmc_clean))
print(gelman.diag(mcmc_clean))

par(mfrow = c(1, 3))
plot(mcmc_clean[, "theta"], main = "Theta", density = FALSE, auto.layout = FALSE)
plot(mcmc_clean[, "sigma"], main = "Sigma", density = FALSE, auto.layout = FALSE)
plot(mcmc_clean[, "alpha"], main = "Alpha", density = FALSE, auto.layout = FALSE)
par(mfrow = c(1, 1))

# Visuals
par(mfrow = c(1, 3))
plot(density(as.matrix(mcmc_clean)[, "theta"]), main="Theta", col="blue", lwd=2)
abline(v=theta_true, col="red", lty=2)

plot(density(as.matrix(mcmc_clean)[, "sigma"]), main="Sigma", col="blue", lwd=2)
abline(v=sigma_true, col="red", lty=2)

plot(density(as.matrix(mcmc_clean)[, "alpha"]), main="Alpha", col="blue", lwd=2)
abline(v=alpha_true, col="red", lty=2)
par(mfrow = c(1, 1))

# Correlation check
mat_all <- as.matrix(mcmc_clean)
cat("\nCorrelations:\n")
print(round(cor(mat_all), 3))

# -----------------------------------------------------------
# 6. RECOVERING G (Finite Sample Diagonal Method)
# -----------------------------------------------------------
mat_all <- as.matrix(mcmc_clean)
y_grid <- seq(0, max(X) * 10, length.out = 200) 
n_sims <- nrow(mat_all)
n_data <- length(X)

gumbel_diag <- function(u, theta, n){
    # METHOD A: Exact
  neg_log_F <- -log(u)
  t_val <- neg_log_F^theta
  t_sum <- n * t_val
  exp(-t_sum^(1/theta))
}

G_exact <- matrix(0, n_sims, length(y_grid))

cat("\nReconstructing G...\n")

for(i in 1:n_sims) {
  th <- mat_all[i, "theta"]
  si <- mat_all[i, "sigma"]
  al <- mat_all[i, "alpha"]
  
  # Mu fixed at 0
  z <- y_grid / si
  
  F_y <- numeric(length(y_grid))
  valid <- z > 0
  F_y[valid] <- pfrechet(z[valid], loc = mu_true, scale = si, shape = al)
  F_y <- pmin(pmax(F_y, 1e-15), 1 - 1e-15)
  
  # METHOD A: Exact
    G[i, ] <- gumbel_diag(F_y, th, n_data)
}

# --- True G ---
z_true <- (y_grid - mu_true) / sigma_true
F_true <- numeric(length(y_grid))
valid_true <- z_true > 0
F_true[valid_true] <- pfrechet(z_true[valid_true], loc = mu_true, scale = sigma_true, shape = alpha_true)

G_true <- gumbel_diag(F_true, theta_true, n_data)

# --- Plotting ---
G_mean <- colMeans(G_exact)
G_lower <- apply(G_exact, 2, quantile, 0.05)
G_upper <- apply(G_exact, 2, quantile, 0.95)

par(mfrow = c(1, 1))
plot(y_grid, G_mean, type = "l", lwd = 2, col = "blue",
     ylim = c(0, 1), 
     main = "Estimated Limit Distribution G",
     xlab = "y", ylab = "P(Mn <= y)")

polygon(c(y_grid, rev(y_grid)), c(G_lower, rev(G_upper)), 
        col = rgb(0, 0, 1, 0.1), border = NA)

lines(y_grid, G_true, col = "red", lwd = 2, lty = 2)