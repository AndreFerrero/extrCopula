# =============================================================================
# Normality check for summary statistics (BSL diagnostic)
# =============================================================================

source("libs/packages.R")

# Load models
source("libs/models/copulas/gumbel.R")
source("libs/models/margins/lognormal.R")
source("libs/models/builders/simulator.R")

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
set.seed(123)

param_map <- list(margin = c("mu", "sigma"), copula = "theta")
simulator <- build_simulator(copula_gumbel, margin_lognormal, param_map)

true_param <- c(mu = 0, sigma = 1, theta = 2)
n_obs <- 1000
n_sim <- 2000   # reasonably large for diagnostics

# Your summary statistics
sum_stats <- function(x) {
  m <- median(x)
  md <- mad(x)
  if (md == 0) md <- 1e-6
  smax <- (max(x) - m) / md
  raw_max <- max(x)
  c(log(m), log(md), log(raw_max))
}

# -----------------------------------------------------------------------------
# Simulate summary statistics
# -----------------------------------------------------------------------------
S <- matrix(NA, n_sim, 3)
colnames(S) <- c("median", "mad", "max")

for (i in 1:n_sim) {
  x <- simulator(true_param, n_obs)
  S[i, ] <- sum_stats(x)
}

pairs(S)
# -----------------------------------------------------------------------------
# Univariate diagnostics
# -----------------------------------------------------------------------------
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))

for (j in 1:3) {
  # Histogram
  hist(S[, j], breaks = 30, probability = TRUE,
       main = paste("Histogram:", colnames(S)[j]),
       xlab = colnames(S)[j])
  curve(
    dnorm(x, mean(S[, j]), sd(S[, j])),
    col = "red", lwd = 2, add = TRUE
  )

  # QQ plot
  qqnorm(S[, j], main = paste("QQ plot:", colnames(S)[j]))
  qqline(S[, j], col = "red", lwd = 2)
}

par(mfrow = c(1, 1))
