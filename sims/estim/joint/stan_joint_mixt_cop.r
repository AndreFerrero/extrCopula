library(copula)
library(rstan)
library(mvtnorm)

# Enable parallel processing for Stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(123)

# -----------------------------------------------------------
# 1. Simulate Data (Same as your Setup)
# -----------------------------------------------------------
n <- 150
theta_true <- 3
mu_true <- 0
sigma_true <- 1

cop <- gumbelCopula(param = theta_true, dim = n)
U_sim <- as.numeric(rCopula(1, cop))
X_sim <- qlnorm(U_sim, meanlog = mu_true, sdlog = sigma_true)
Y_sim <- log(X_sim) # Working on Log-Scale

# -----------------------------------------------------------
# 2. Prepare for Stan
# -----------------------------------------------------------
stan_data <- list(
  D = n,
  Y = Y_sim
)

# -----------------------------------------------------------
# 3. Compile and Run
# -----------------------------------------------------------

joint_dir <- 
mod <- stan_model(here(joint_dir, "bayes_mixture_copula.stan"))

# Run NUTS Sampler
# NUTS is much more efficient than Metropolis-Hastings for this
fit <- sampling(mod, 
                data = stan_data, 
                chains = 4, 
                iter = 2000, 
                warmup = 1000,
                control = list(adapt_delta = 0.95)) # Higher delta for safety

res_dir <- here(joint_dir, "res")
save(fit, file = here(res_dir, "stan_fit.Rdata"))

# -----------------------------------------------------------
# 4. Diagnostics & Results
# -----------------------------------------------------------
print(fit, pars = c("theta", "w", "mu", "sigma"))

# Traceplot
traceplot(fit, pars = c("theta"))

# Density Plot
stan_dens(fit, pars = "theta", separate_chains = TRUE) + 
  geom_vline(xintercept = theta_true, color = "red", linetype = "dashed")

# Extract samples
post_samples <- rstan::extract(fit)
cat("\nMean Posterior Theta:", mean(post_samples$theta), "\n")
cat("80% CI:", quantile(post_samples$theta, probs = c(0.1, 0.9)), "\n")
