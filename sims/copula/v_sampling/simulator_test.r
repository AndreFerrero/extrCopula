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
  copula = copula_clayton,
  margin = margin_lognormal,
  param_map = param_map
)

# ------------------------------
# 5. Test simulation
# ------------------------------
param_test <- c(mu = 0, sigma = 1, theta = 3)
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

hist(X_sim, breaks = 30,
     main = "Simulated Gumbel + LogNormal",
     xlab = "X",
     col = "skyblue")
