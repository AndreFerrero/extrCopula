# =============================================================================
# Test script: small MCMC run with Gumbel copula and Lognormal margin
# =============================================================================

# Load packages
source("libs/packages.R")

# Load models
source("libs/models/copulas/gumbel.R")
source("libs/models/margins/lognormal.R")

# Load builders
source("libs/models/builders/simulator.R")
source("libs/models/builders/logposterior.R")
source("libs/models/builders/semibsl_logposterior.R")


# Load MCMC machinery
source("libs/mcmc/run_chain.R")
source("libs/mcmc/engines/metropolis_hastings.R")
source("libs/mcmc/proposals/gaussian_rw.R")
source("libs/mcmc/adaptation/none.R")
source("libs/mcmc/adaptation/haario.R")
source("libs/mcmc/adaptation/robbins_monro.R")
# =============================================================================
# 1. Simulate data using the simulator
# =============================================================================
param_map <- list(margin = c("mu", "sigma"), copula = "theta")
simulator <- build_simulator(copula_gumbel, margin_lognormal, param_map)

set.seed(123)
true_param <- c(mu = 0, sigma = 1, theta = 3)
n_obs <- 1000
X <- simulator(true_param, n_obs)

# =============================================================================
# 2. Build log-posterior
# =============================================================================
to_unconstrained <- function(param) {
  c(param["mu"], log(param["sigma"]), log(param["theta"] - 1))
}

from_unconstrained <- function(phi) {
  param <- c(
    mu    = phi[1],
    sigma = exp(phi[2]),
    theta = exp(phi[3]) + 1
  )
  # Ensure names are exactly what we want
  names(param) <- c("mu", "sigma", "theta")
  param
}


log_jacobian <- function(phi) phi[2] + phi[3]

logpost <- build_semibsl_logposterior(
  simulator = fn_sim_gumbel,
  sum_stats = fn_sum_stats,
  n_sim = 200,
  log_prior = fn_log_prior,
  transform = list(
    to_unconstrained = to_unconstrained,
    to_natural = to_natural,
    log_jacobian = log_jacobian
  ),
  obs_data = X
)

# =============================================================================
# 3. Define proposal and run chain
# =============================================================================
phi_init <- to_unconstrained(c(mu = 0.5, sigma = 0.8, theta = 2.5))
proposal <- proposal_gaussian_rw(Sigma = diag(0.01, 3))

res <- run_chain(
  log_target = logpost,
  init = phi_init,
  n_iter = 20000,
  proposal = proposal,
  adapt = adapt_robbins_monro(target_accept = 0.234)
)

# Convert samples back to natural space
samples_natural <- t(apply(res$samples, 1, from_unconstrained))

# =============================================================================
# 4. Quick summaries
# =============================================================================
cat("Acceptance rate:", res$accept_rate, "\n")
cat("Posterior means:\n")
print(colMeans(samples_natural[10000:20000, ]))

# Traceplot

# Arrange 1 row and 3 columns
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))  # smaller margins for tighter plots

plot(samples_natural[, "mu"], type = "l",
     ylab = "mu", xlab = "Iteration", main = "Trace of mu")

plot(samples_natural[, "sigma"], type = "l",
     ylab = "sigma", xlab = "Iteration", main = "Trace of sigma")

plot(samples_natural[, "theta"], type = "l",
     ylab = "theta", xlab = "Iteration", main = "Trace of theta")

# Reset par
par(mfrow = c(1,1))
