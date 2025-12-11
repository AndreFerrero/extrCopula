# ==============================================================================
# Semi-Parametric Bayesian Synthetic Likelihood with Adaptive Block Metropolis
# ==============================================================================

library(parallel)
library(stabledist)
library(mvtnorm)
library(MASS)
library(coda)

# ==============================================================================
# 2. SemiBSL Functions
# ==============================================================================

# 2.1 Gumbel Copula Simulation (V-stable)
rGumbV <- function(n, theta) {
  val_gamma <- (cos(pi / (2 * theta)))^theta
  V <- rstable(n = 1, alpha = 1 / theta, beta = 1, gamma = val_gamma, delta = 0, pm = 1)
  if (V <= 0) V <- .Machine$double.eps
  E <- rexp(n)
  exp(-(E / V)^(1 / theta))
}

# 2.2 Simulator
fn_sim_gumbel <- function(param_vec, n_obs) {
  mu <- param_vec[1]
  sigma <- param_vec[2]
  theta <- param_vec[3] # Copula parameter
  if (sigma <= 0 || theta < 1) {
    return(NA)
  }

  U <- tryCatch(rGumbV(n_obs, theta), error = function(e) NULL)
  if (is.null(U) || any(is.na(U))) {
    return(NA)
  }

  X <- qlnorm(U, meanlog = mu, sdlog = sigma)
  if (any(is.infinite(X))) {
    return(NA)
  }
  return(X)
}

# 2.3 Summary statistics
fn_sum_stats <- function(x) {
  if (any(is.na(x))) {
    return(rep(NA, 5))
  }
  m <- median(x)
  md <- mad(x)
  q_diff <- diff(quantile(x, probs = c(0.05, 0.95)))
  tail_ratio <- (quantile(x, 0.95) - m) / (m - quantile(x, 0.05))
  if (md == 0) md <- 1e-6
  s_max <- (max(x) - m) / md
  c(m, md, q_diff, tail_ratio, s_max)
}

# 2.4 Prior
fn_log_prior <- function(param_vec) {
  mu <- param_vec[1]
  sigma <- param_vec[2]
  theta <- param_vec[3]

  if (sigma <= 0 || theta <= 1) {
    return(-Inf)
  }

  lp_mu <- dnorm(mu, mean = 0, sd = 10, log = TRUE)
  lp_sigma <- dcauchy(sigma, location = 0, scale = 2.5, log = TRUE)
  lp_theta <- dgamma(theta - 1, shape = 2, rate = 0.5, log = TRUE)

  lp_mu + lp_sigma + lp_theta
}

# 2.5 KDE for semi-parametric marginal likelihood
get_kde_estimates <- function(sim_data, obs_val) {
  n <- length(sim_data)
  sd_sim <- sd(sim_data)
  iqr_sim <- IQR(sim_data)
  h <- 0.9 * min(sd_sim, iqr_sim / 1.34) * n^(-0.2)
  if (h == 0) h <- 1e-6
  z_scores <- (obs_val - sim_data) / h
  log_vals <- dnorm(z_scores, log = TRUE)
  max_val <- max(log_vals)
  log_pdf <- max_val + log(sum(exp(log_vals - max_val))) - log(n) - log(h)
  u_obs <- mean(pnorm(z_scores))
  u_sim <- rank(sim_data) / (n + 1)
  list(log_pdf = log_pdf, u_obs = u_obs, u_sim = u_sim)
}

# 2.6 SemiBSL Log-Likelihood
semi_bsl_loglik <- function(param_vec, y_sum, n_sim, n_obs) {
  d <- length(y_sum)
  sim_stats <- matrix(NA, nrow = n_sim, ncol = d)

  for (i in 1:n_sim) {
    x_sim <- fn_sim_gumbel(param_vec, n_obs)
    if (any(is.na(x_sim))) {
      return(-Inf)
    }
    sim_stats[i, ] <- fn_sum_stats(x_sim)
  }

  total_loglik <- 0
  eta_obs <- numeric(d)
  eta_sim <- matrix(NA, nrow = n_sim, ncol = d)

  for (j in 1:d) {
    kde_res <- get_kde_estimates(sim_stats[, j], y_sum[j])
    total_loglik <- total_loglik + kde_res$log_pdf
    eta_obs[j] <- qnorm(max(min(kde_res$u_obs, 1 - 1e-6), 1e-6))
    eta_sim[, j] <- qnorm(kde_res$u_sim)
  }

  R_hat <- cor(eta_sim)
  if (det(R_hat) <= 1e-9) R_hat <- R_hat + diag(1e-6, d)
  R_inv <- tryCatch(solve(R_hat), error = function(e) MASS::ginv(R_hat))
  log_det_R <- as.numeric(determinant(R_hat, logarithm = TRUE)$modulus)
  quadratic_term <- as.numeric(t(eta_obs) %*% (R_inv - diag(1, d)) %*% eta_obs)
  copula_loglik <- -0.5 * log_det_R - 0.5 * quadratic_term

  total_loglik + copula_loglik
}

# 2.7 Transformations
to_unconstrained <- function(param_vec) c(param_vec[1], log(param_vec[2]), log(param_vec[3] - 1))
to_natural <- function(phi) c(phi[1], exp(phi[2]), exp(phi[3]) + 1)
log_jacobian <- function(phi) phi[2] + phi[3]

# ==============================================================================
# 3. Adaptive Block Metropolis for SemiBSL
# ==============================================================================

run_worker_bsl <- function(y_obs_sum, n_mcmc, burn_in, n_sim, n_obs, init_param, chain_id = 1) {
  d_param <- 3
  chain <- matrix(NA, nrow = n_mcmc, ncol = d_param)
  colnames(chain) <- c("mu", "sigma", "theta")
  acc <- 0

  curr_phi <- to_unconstrained(init_param)
  curr_param <- init_param
  curr_ll <- semi_bsl_loglik(curr_param, y_obs_sum, n_sim, n_obs)
  curr_lp <- fn_log_prior(curr_param)
  curr_post <- curr_ll + curr_lp + log_jacobian(curr_phi)
  if (!is.finite(curr_post)) stop("Initial param invalid")

  cov_mat <- diag(0.01, d_param)
  sd_scale <- (2.38^2) / d_param
  eps <- 1e-6 * diag(d_param)

  pb <- txtProgressBar(min = 0, max = n_mcmc, style = 3)

  for (i in 1:n_mcmc) {
    prop_phi <- as.vector(rmvnorm(1, mean = curr_phi, sigma = cov_mat))
    prop_param <- to_natural(prop_phi)

    prop_lp <- fn_log_prior(prop_param)
    prop_jac <- log_jacobian(prop_phi)

    if (prop_lp == -Inf) {
      log_alpha <- -Inf
    } else {
      prop_ll <- semi_bsl_loglik(prop_param, y_obs_sum, n_sim, n_obs)
      prop_post <- prop_ll + prop_lp + prop_jac
      log_alpha <- prop_post - curr_post
    }

    if (log(runif(1)) < log_alpha) {
      curr_phi <- prop_phi
      curr_param <- prop_param
      curr_post <- prop_post
      curr_ll <- prop_lp
      acc <- acc + 1
    }

    chain[i, ] <- curr_param

    # Adaptive covariance during burn-in
    if (i > 50 && i <= burn_in && i %% 50 == 0) {
      hist_phi <- t(apply(chain[1:i, ], 1, to_unconstrained))
      cov_mat <- sd_scale * cov(hist_phi) + eps
    }

    # Progress
    setTxtProgressBar(pb, i)
    if (i %% 50 == 0) {
      cat(sprintf("Chain %d | Iter %d/%d | Acc: %.2f%%\n", chain_id, i, n_mcmc, 100 * acc / i))
    }
  }
  close(pb)
  cat(sprintf("Chain %d | Final Acceptance Rate: %.2f%%\n", chain_id, 100 * acc / n_mcmc))
  mcmc(chain)
}

# ==============================================================================
# 4. Run Parallel Chains
# ==============================================================================
set.seed(123)

# ==============================================================================
# 1. Simulate Data (Log-Normal + Gumbel Copula)
# ==============================================================================
set.seed(123)
n <- 1000
mu_true <- 0
sigma_true <- 1
theta_true <- 3

param_true <- c(mu_true, sigma_true, theta_true)

X <- fn_sim_gumbel(param_true, n)

y_sum <- fn_sum_stats(X)

n_chains <- 4
n_mcmc <- 10000
burn_in <- floor(n_mcmc / 2)
n_sim_bsl <- 800

inits_list <- list(
  c(0.5, 0.8, 1.8),
  c(-0.5, 1.2, 2.5),
  c(0.1, 1.0, 1.8),
  c(-0.1, 1.5, 3.0)
)

cl <- makeCluster(n_chains, outfile = "")
clusterEvalQ(cl, {
  library(mvtnorm)
  library(MASS)
  library(stabledist)
})

clusterExport(cl, varlist = c(
  "X", "y_sum", "n_mcmc", "burn_in", "n_sim_bsl", "inits_list",
  "fn_sum_stats", "fn_sim_gumbel", "rGumbV",
  "fn_log_prior", "semi_bsl_loglik",
  "to_unconstrained", "to_natural", "log_jacobian",
  "run_worker_bsl", "get_kde_estimates"
))

cat("Starting Parallel SemiBSL MCMC...\n")
results <- parLapply(cl, 1:n_chains, function(cid) {
  run_worker_bsl(
    y_obs_sum = y_sum,
    n_mcmc = n_mcmc,
    burn_in = burn_in,
    n_sim = n_sim_bsl,
    n_obs = length(X),
    init_param = inits_list[[cid]],
    chain_id = cid
  )
})
stopCluster(cl)

sbi_dir <- here("sims", "estim", "joint", "SBI")
res_dir <- here(sbi_dir, "res")
save(results, file = here(res_dir, "adaptstop_semiBSL_lognorm.Rdata"))

# ==============================================================================
# 5. Diagnostics
# ==============================================================================
load(here(res_dir, "adaptstop_semiBSL_lognorm.Rdata"))
mcmc_list <- mcmc.list(results)
mcmc_clean <- window(mcmc_list, start = burn_in + 1)

summary(mcmc_clean)
gelman.diag(mcmc_clean)
effectiveSize(mcmc_clean)

# Traceplots
par(mfrow = c(1, 3))
plot(mcmc_clean[, "theta"], main = "Copula Theta")
plot(mcmc_clean[, "sigma"], main = "Sigma")
plot(mcmc_clean[, "mu"], main = "Mu")
par(mfrow = c(1, 1))
