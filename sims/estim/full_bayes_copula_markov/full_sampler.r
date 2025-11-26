# Load required packages
library(copula)
library(MASS)   # for normal distributions in mixture
library(rstan)  # for GEV posterior sampling

exp_folder <- here("sims", "estim", "full_bayes_copula_markov")

# -----------------------------
# 1. Model parameters
# -----------------------------
n <- 500           # number of observations
M <- 1000          # number of posterior iterations
K_max <- 20        # max clusters for DP approximation
alpha0 <- 1        # DP concentration parameter
G0_mu <- 0         # base measure for mixture (e.g., normal)
G0_sigma <- 1

# Hyperpriors for copula
alpha_prior_a <- 1
alpha_prior_b <- 1

library(copula)

set.seed(123)

# Parameters
n <- 500       # length of chain
burn <- 100    # burn-in
alpha <- 2.5   # Gumbel copula parameter
gumb <- gumbelCopula(alpha)

# Initialize
U <- numeric(n + burn)
U[1] <- runif(1)  # start uniform

# Simulate Markov chain
for (t in 1:(n + burn - 1)) {
  V <- runif(1)
  U[t + 1] <- cCopula(cbind(U[t], V), copula = gumb, inverse = TRUE)[2]
}

# Drop burn-in
U <- U[(burn + 1):(n + burn)]

# Transform to desired margins, e.g., standard normal
X <- qnorm(U)  # or qfrechet(U) if you want heavy tails

# -----------------------------
# 2. Initialize storage
# -----------------------------
F_posterior <- vector("list", M)           # stores DP mixture parameters
U_chains <- matrix(NA, nrow=n, ncol=M)     # pseudo-observations
alpha_chain <- numeric(M)
theta_chain <- numeric(M)

# -----------------------------
# 3. Stage 1: Sample F, alpha, theta
# -----------------------------
for (m in 1:M) {
  
  # ----- A. Neal's Algorithm 8 for DP mixture F -----
  # Placeholder: assume normal mixture; assign cluster labels randomly
  # For a real implementation, follow Neal's Algorithm 8:
  # 1. sample cluster assignments for each X_t
  # 2. sample cluster parameters
  # Here we just illustrate the structure
  cluster_assign <- sample(1:K_max, n, replace=TRUE)
  cluster_params <- mapply(function(k) rnorm(1, G0_mu, G0_sigma), 1:K_max)
  
  F_posterior[[m]] <- list(clusters=cluster_assign, params=cluster_params)
  
  # ----- B. Compute pseudo-observations U_t -----
  # For normal mixture, compute cdf of mixture for each X_t
  U_t <- numeric(n)
  for (i in 1:n) {
    k <- cluster_assign[i]
    mu_k <- cluster_params[k]
    sigma_k <- 1  # fixed for simplicity
    U_t[i] <- pnorm(X[i], mean=mu_k, sd=sigma_k)  # mixture cdf approx
  }
  U_chains[, m] <- U_t
  
  # ----- C. Sample copula parameter alpha (Gumbel) -----
  # Compute log-likelihood of Gumbel copula
  loglik <- function(alpha) {
    cop <- gumbelCopula(alpha)
    sum(dCopula(cbind(U_t[-n], U_t[-1]), copula=cop, log=TRUE))
  }
  # MH sampler for alpha
  alpha_curr <- ifelse(m==1, 2, alpha_chain[m-1])
  alpha_prop <- rnorm(1, alpha_curr, 0.1)
  if (alpha_prop > 1) {  # Gumbel alpha > 1
    log_accept <- loglik(alpha_prop) - loglik(alpha_curr)
    if (log(runif(1)) < log_accept) alpha_curr <- alpha_prop
  }
  alpha_chain[m] <- alpha_curr
  
  # ----- D. Compute extremal index theta -----
  # Analytic formula for Gumbel copula: theta = 2^(1/alpha) - 1
  theta_chain[m] <- 2^(1/alpha_curr) - 1
  
  cat("Iteration", m, "done\n")
}

# Compute empirical CDF
ecdf_X <- ecdf(X)

plot(ecdf_X, main="Posterior draws of F vs empirical CDF", xlab="X", ylab="CDF")

# Overlay a few posterior draws
for (m in sample(1:M, 50)) {  # take 50 draws at random
  lines(sort(X), sort(U_chains[, m]), col=rgb(0,0,1,0.2))
}

# -----------------------------
# 4. Stage 2: Fit GEV to block maxima
# -----------------------------
# Define block size
block_size <- 50
n_blocks <- floor(n / block_size)
M_blocks <- sapply(1:n_blocks, function(b) {
  max(X[((b-1)*block_size + 1):(b*block_size)])
})

# Stan model for GEV
stan_data <- list(N=n_blocks, M=M_blocks)

stan_file <- here(exp_folder, "gev.stan")
sm <- stan_model(stan_file)

fit_gev <- sampling(sm,
                data = stan_data,
                iter = 3000,
                warmup = 1000,
                chains = 4,
                control = list(adapt_delta = 0.99, max_treedepth = 12),
                seed = 2025)
gev_post <- extract(fit_gev)

# -----------------------------
# 5. Combine theta with distorted GEV
# -----------------------------
n_post <- length(gev_post$mu)
mu_iid <- numeric(n_post)
sigma_iid <- numeric(n_post)
xi_iid <- numeric(n_post)

for (i in 1:n_post) {
  theta_i <- theta_chain[sample(1:M, 1)]  # randomly pair posterior draws
  mu_d <- gev_post$mu[i]
  sigma_d <- gev_post$sigma[i]
  xi_d <- gev_post$xi[i]
  
  mu_iid[i] <- mu_d - sigma_d / xi_d * (1 - theta_i^xi_d)
  sigma_iid[i] <- sigma_d * theta_i^xi_d
  xi_iid[i] <- xi_d
}

# Result: posterior chains of i.i.d. GEV parameters
