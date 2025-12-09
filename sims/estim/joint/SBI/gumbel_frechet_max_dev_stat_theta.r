library(mvtnorm)
library(compiler)

set.seed(123)

# -----------------------------------------------------------
# 1. FAST SIMULATOR (Copied from previous steps)
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
  X <- mu + sigma * (-log(U_vec))^(-1/alpha) # Frechet Margin
  return(X)
})

# -----------------------------------------------------------
# 2. DEFINING THE STATISTIC
# -----------------------------------------------------------
calc_dep_stat <- function(x) {
  
  # The "Normalized Tail Reach"
  stat <- (max(x) - median(x)) / mad(x)
  return(stat)
}

# -----------------------------------------------------------
# 3. THE EXPERIMENT
# -----------------------------------------------------------
# Fixed Margins
mu <- 10; sigma <- 5; alpha <- 3
n_dim <- 1000

# Grid of Thetas to test
theta_grid <- c(1.1, 1.5, 2, 3, 5)
n_reps <- 200 # Simulations per theta

results <- list()

for(th in theta_grid) {
  stats <- numeric(n_reps)
  for(i in 1:n_reps) {
    x <- sim_gumbel_fast(th, mu, sigma, alpha, n=n_dim)
    stats[i] <- calc_dep_stat(x)
  }
  results[[paste0("theta_", th)]] <- stats
}

# -----------------------------------------------------------
# 4. VISUALIZATION
# -----------------------------------------------------------
boxplot(results, 
        main="Sensitivity of 'Max-Median / MAD' to Theta",
        xlab="Theta (Dependence)",
        ylab="Statistic Value",
        col="lightblue",
        las=1) # Rotate labels

# Add a smoothing line through the medians
medians <- sapply(results, median)
lines(1:length(theta_grid), medians, col="red", lwd=2, type="b")

# Check monotonicity
print(medians)
