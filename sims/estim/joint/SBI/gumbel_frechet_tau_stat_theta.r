library(mvtnorm)
library(compiler)

set.seed(123)

# -----------------------------------------------------------
# 1. GENERATOR (Your fast implementation)
# -----------------------------------------------------------
sim_gumbel_fast <- cmpfun(function(theta, mu, sigma, alpha, n=150) {
  a <- 1/theta
  U_pi <- runif(1, -pi/2, pi/2)
  W    <- rexp(1)
  t1 <- sin(a * (U_pi + pi/2))
  t2 <- cos(U_pi)^(1/a)
  t3 <- (cos(U_pi - a * (U_pi + pi/2)) / W)^((1-a)/a)
  V <- (t1 / t2) * t3 * (cos(pi * a / 2))^(theta)
  
  E <- rexp(n)
  U_vec <- exp( - (E / V)^a )
  X <- mu + sigma * (-log(U_vec))^(-1/alpha) 
  return(X)
})

# -----------------------------------------------------------
# 2. DEFINING THE STATISTIC (Kendall's Tau)
# -----------------------------------------------------------
calc_tau_stat <- function(x) {
  # We estimate pairwise correlation using lag-1 pairs
  # (X_1, X_2), (X_2, X_3)...
  # Since the copula is exchangeable, any pair works.
  
  val <- cor(x[-length(x)], x[-1], method="kendall")
  return(val)
}

# -----------------------------------------------------------
# 3. THE EXPERIMENT
# -----------------------------------------------------------
mu <- 10; sigma <- 5; alpha <- 3
n_dim <- 1000
theta_grid <- c(1.1, 1.5, 2, 3, 5, 8)
n_reps <- 200

results_tau <- list()

for(th in theta_grid) {
  stats <- numeric(n_reps)
  for(i in 1:n_reps) {
    x <- sim_gumbel_fast(th, mu, sigma, alpha, n=n_dim)
    stats[i] <- calc_tau_stat(x)
  }
  results_tau[[paste0("theta_", th)]] <- stats
}

# -----------------------------------------------------------
# 4. VISUALIZATION
# -----------------------------------------------------------
# Theoretical Tau for Gumbel is 1 - 1/theta
theoretical_tau <- 1 - 1/theta_grid

boxplot(results_tau, 
        main="Sensitivity of Kendall's Tau to Theta",
        xlab="Theta (Dependence)",
        ylab="Kendall's Tau",
        col="#ffcc99",
        las=1)

# Add trend line
medians <- sapply(results_tau, median)
lines(1:length(theta_grid), medians, col="darkred", lwd=2, type="b")

# Add Theoretical Truth (Green) to verify calibration
points(1:length(theta_grid), theoretical_tau, col="darkgreen", pch=4, lwd=3, cex=1.5)
legend("topleft", legend=c("Simulation Median", "Theory (1 - 1/theta)"), 
       col=c("darkred", "darkgreen"), lty=c(1, NA), pch=c(1, 4), lwd=c(2,2))