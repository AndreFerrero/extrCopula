# ------------------------------
# 1. Load required packages
# ------------------------------
source("libs/packages.R")

# ------------------------------
# 2. Source your libs
# ------------------------------
source("libs/models/copulas/gumbel.R")     # defines copula_gumbel
source("libs/models/copulas/clayton.r")    # defines copula_clayton
source("libs/models/margins/lognormal.R")  # defines margin_lognormal
source("libs/models/builders/simulator.R") # defines build_simulator()

# ------------------------------
# 3. Define parameter mapping
# ------------------------------
# param vector will be named, e.g., c(mu, sigma, theta)
param_map <- list(
  margin = c("mu", "sigma"),
  copula = "theta"
)

# ------------------------------
# 4. Build the simulator
# ------------------------------
simulator <- build_simulator(
  copula = copula_gumbel,
  margin = margin_lognormal,
  param_map = param_map
)

# ------------------------------
# 5. Test simulation
# ------------------------------
param_test <- c(mu = 0, sigma = 1, theta = 1)
n_test <- 1000

X_sim <- simulator(param_test, n_test)

# ------------------------------
# 6. Quick sanity checks
# ------------------------------
if (all(is.na(X_sim))) {
  stop("Simulator returned NA: check copula parameters or numerical stability")
}

cat("First 10 simulated values:\n")
print(head(X_sim, 10))

cat("\nSummary statistics:\n")
print(summary(X_sim))

# ------------------------------
# Histogram of simulated data
# ------------------------------
hist(X_sim, breaks = 30,
     main = "Simulated Gumbel + LogNormal",
     xlab = "X",
     col = "skyblue",
     probability = TRUE)  # scale histogram to density

# ------------------------------
# Overlay lognormal density
# ------------------------------
x_vals <- seq(min(X_sim), max(X_sim), length.out = 200)
lines(x_vals, dlnorm(x_vals, meanlog = param_test["mu"], sdlog = param_test["sigma"]),
      col = "red", lwd = 2, lty = 2)  # dashed red line

# ------------------------------
# Overlay KDE of simulated data
# ------------------------------
dens <- density(X_sim)
lines(dens, col = "darkgreen", lwd = 2)  # smooth green line

# ------------------------------
# Add legend
# ------------------------------
legend("topright", legend = c("Histogram", "LogNormal density", "KDE"),
       col = c("skyblue", "red", "darkgreen"), lwd = c(10, 2, 2),
       lty = c(1, 2, 1), bty = "n")

theta_vec <- c(1, 1.5, 2, 3, 5)
n_sim <- 1000
results <- list()

for(theta in theta_vec) {
  param <- c(mu = 0, sigma = 1, theta = theta)
  X <- simulator(param, n_sim)
  
  # Store results
  results[[as.character(theta)]] <- list(
    data = X,
    mean = mean(X),
    sd = sd(X),
    quantiles = quantile(X, probs = c(0.05, 0.5, 0.95)),
    # correlation of components if X is multivariate
    cor = if(is.matrix(X)) cor(X) else NA
  )
}

# Example: plot KDEs over theta
plot(density(results[["1"]]$data), col = "black", lwd = 2,
     main = "Effect of Theta on Distribution",
     xlab = "X", ylim = c(0, 0.5))
cols <- c("red","blue","green","purple")
for(i in 2:length(theta_vec)){
  lines(density(results[[as.character(theta_vec[i])]]$data), 
        col = cols[i-1], lwd = 2)
}
legend("topright", legend = paste0("theta=", theta_vec),
       col = c("black", cols), lwd = 2)
