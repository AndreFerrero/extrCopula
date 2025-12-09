library(copula)
library(coda)
library(mvtnorm)
library(future)
library(future.apply)

set.seed(123)

# -----------------------------------------------------------
# 1. Simulate Data (Log-Normal Margins)
# -----------------------------------------------------------
n <- 100000
theta_true <- 3

# Log-Normal Parameters (corresponding to Normal on log-scale)
mu_true <- 0
sigma_true <- 1

paste("1. Simulating Data (Log-Normal Margins, n =", n, ")...\n")

# A. Dependence
cop <- gumbelCopula(param = theta_true, dim = n)
U <- as.numeric(rCopula(1, cop))

# B. Margins (Log-Normal)
# X are the raw observations
X <- qlnorm(U, meanlog = mu_true, sdlog = sigma_true)

# -----------------------------------------------------------
# 2. Model Definitions (Log-Normal + Gumbel)
# -----------------------------------------------------------

log_posterior <- function(param_vec, data) {
    # Unpack parameters (Only 3 now: theta, mu, sigma)
    theta <- param_vec[1]
    mu <- param_vec[2]
    sigma <- param_vec[3]

    # --- A. PRIORS ---

    # Copula: Gamma(2,1) shifted
    if (theta <= 1.01 || theta > 20) {
        return(-Inf)
    }
    lp_theta <- dgamma(theta - 1, shape = 2, rate = 1, log = TRUE)

    # Location (mu): Vague Normal
    # Note: Ideally centered near mean(log(data))
    lp_mu <- dnorm(mu, mean = 0, sd = 10, log = TRUE)

    # Scale (sigma): Boundary Avoiding
    # Even though LogNormal is identifiable, keeping the lower bound
    # is good practice in N=1 joint estimation.
    if (sigma <= 0.01) {
        return(-Inf)
    }
    lp_sigma <- dlnorm(sigma, meanlog = 0, sdlog = 1, log = TRUE)

    # --- B. MARGINAL LIKELIHOOD (Log-Normal) ---
    # dlnorm calculates density of Log-Normal
    dens_vals <- dlnorm(data, meanlog = mu, sdlog = sigma, log = TRUE)
    ll_margin <- sum(dens_vals)

    # --- C. COPULA LIKELIHOOD (Full Density) ---
    # 1. Transform Data to U using plnorm (Log-Normal CDF)
    u_hat <- plnorm(data, meanlog = mu, sdlog = sigma)

    # 2. Clamp for numerical safety
    u_hat <- pmin(pmax(u_hat, 1e-8), 1 - 1e-8)

    # 3. Calculate Full Likelihood
    ld_cop <- dCopula(u_hat,
        copula = gumbelCopula(theta, dim = length(data)),
        log = TRUE
    )

    # 4. Check for numerical failure
    if (is.na(ld_cop) || !is.finite(ld_cop)) {
        return(-Inf)
    }

    # Total Posterior
    return(lp_theta + lp_mu + lp_sigma + ll_margin + ld_cop)
}

# -----------------------------------------------------------
# 3. Worker: Block Adaptive Metropolis
# -----------------------------------------------------------
run_worker <- function(data, n_iter, init_vals, chain_id, progress_dir) {
    library(copula)
    library(mvtnorm)
    library(coda)

    p_names <- c("theta", "mu", "sigma")
    n_par <- length(p_names)
    chain <- matrix(NA, nrow = n_iter, ncol = n_par)
    colnames(chain) <- p_names

    curr_par <- init_vals

    # Robust Start
    curr_lp <- log_posterior(curr_par, data)
    if (!is.finite(curr_lp)) {
        stop("Worker failed: Initial values yield -Inf posterior.")
    }

    # --- ADAPTIVE METROPOLIS SETTINGS ---
    cov_mat <- diag(rep(0.01, n_par))
    # Optimal Scaling for d=3
    sd_scale <- (2.38^2) / n_par
    eps <- 1e-6 * diag(n_par)

    adapt_start <- 500
    total_accepts <- 0
    status_file <- file.path(progress_dir, paste0("status_", chain_id, ".txt"))

    for (i in 1:n_iter) {
        # 1. BLOCK PROPOSAL
        prop_val <- rmvnorm(1, mean = curr_par, sigma = cov_mat)[1, ]

        prop_lp <- log_posterior(prop_val, data)
        ratio <- prop_lp - curr_lp

        # 2. ACCEPT/REJECT
        if (is.finite(ratio) && log(runif(1)) < ratio) {
            curr_par <- prop_val
            curr_lp <- prop_lp
            total_accepts <- total_accepts + 1
        }
        chain[i, ] <- curr_par

        # 3. ADAPTATION
        if (i > adapt_start && i %% 50 == 0) {
            emp_cov <- cov(chain[1:i, ])
            cov_mat <- (sd_scale * emp_cov) + eps
        }

        # 4. DASHBOARD
        if (i %% 50 == 0 || i == n_iter) {
            rate <- (total_accepts / i) * 100
            pct <- as.integer((i / n_iter) * 100)
            try(
                {
                    writeLines(sprintf("%d|%.1f", pct, rate), status_file)
                },
                silent = TRUE
            )
        }
    }
    mcmc(chain)
}

# -----------------------------------------------------------
# 4. Execution
# -----------------------------------------------------------
n_chains <- 3
plan(multisession, workers = n_chains)
progress_dir <- tempdir()
file.remove(list.files(progress_dir, pattern = "status_", full.names = TRUE))

n_iter <- 10000

inits_list <- list()
for (c in 1:n_chains) {
    inits_list[[c]] <- c(
        theta = runif(1, 1.5, 5),
        mu = mean(log(X)), # Initialize near sample moments of log data
        sigma = sd(log(X))
    )
}

futures_list <- list()
for (c in 1:n_chains) {
    futures_list[[c]] <- future(
        {
            run_worker(X, n_iter, inits_list[[c]], c, progress_dir)
        },
        seed = TRUE
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

# -----------------------------------------------------------
# 5. Analysis
# -----------------------------------------------------------
res_dir <- here("sims", "estim", "joint", "res")

# save(futures_list, file = here(res_dir, "lognormal_bayes_chains.Rdata"))
load(here(res_dir, "lognormal_bayes_chains.Rdata"))

burn_in <- n_iter / 2

chains_list <- value(futures_list)
mcmc_obj <- mcmc.list(chains_list)
mcmc_clean <- window(mcmc_obj, start = burn_in + 1, thin = 5)

summary(mcmc_clean)


# DIAGNOSTICS
print(gelman.diag(mcmc_clean))
print(effectiveSize(mcmc_clean))

# TRACEPLOTS
par(mfrow = c(1, 3))
plot(mcmc_clean[, "theta"], main = "Theta", density = FALSE, auto.layout = FALSE)
plot(mcmc_clean[, "sigma"], main = "Sigma", density = FALSE, auto.layout = FALSE)
plot(mcmc_clean[, "mu"], main = "Mu", density = FALSE, auto.layout = FALSE)
par(mfrow = c(1, 1))

# ESTIMATES
paste("True Theta:", theta_true)
paste("Posterior Mean Theta:", round(mean(as.matrix(mcmc_clean[, "theta"])), 3))
paste("Posterior Median Theta:", round(median(as.matrix(mcmc_clean[, "theta"])), 3))
paste("90% CI:", round(quantile(as.matrix(mcmc_clean[, "theta"]), probs = c(0.1, 0.9)), 3))

paste("True Mu:", mu_true)
paste("Posterior Mean Mu:", round(mean(as.matrix(mcmc_clean[, "mu"])), 3))
paste("Posterior Median Mu:", round(median(as.matrix(mcmc_clean[, "mu"])), 3))

paste("True Sigma:", sigma_true)
paste("Posterior Mean Sigma:", round(mean(as.matrix(mcmc_clean[, "sigma"])), 3))
paste("Posterior Median Sigma:", round(median(as.matrix(mcmc_clean[, "sigma"])), 3))

####################
# DENSITIES
par(mfrow = c(1, 3))
theta_dens_est <- density(as.matrix(mcmc_clean)[, "theta"])
plot(theta_dens_est,
    main = "Posterior Density of Theta",
    lwd = 2, col = "blue", xlab = expression(theta)
)
abline(v = theta_true, col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)

mu_dens_est <- density(as.matrix(mcmc_clean)[, "mu"])
plot(mu_dens_est,
    main = "Posterior Density of Mu",
    lwd = 2, col = "blue", xlab = expression(mu)
)
abline(v = mu_true, col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)

sigma_dens_est <- density(as.matrix(mcmc_clean)[, "sigma"])
plot(sigma_dens_est,
    main = "Posterior Density of Sigma",
    lwd = 2, col = "blue", xlab = expression(sigma)
)
abline(v = sigma_true, col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)
par(mfrow = c(1, 1))


# 1. Define the range you want to view
# (Make sure these numbers are within the range of your current mcmc_clean)
start_iter <- 1
end_iter <- 1000

# 2. Create a subset object using window()
mcmc_zoom <- window(mcmc_obj, start = start_iter, end = end_iter)

# 3. Plot
par(mfrow = c(1, 3))

plot(mcmc_zoom[, "theta"],
    main = "Theta",
    density = FALSE, smooth = FALSE, auto.layout = FALSE
)
plot(mcmc_zoom[, "sigma"],
    main = "Sigma",
    density = FALSE, smooth = FALSE, auto.layout = FALSE
)
plot(mcmc_zoom[, "mu"],
    main = "Mu",
    density = FALSE, smooth = FALSE, auto.layout = FALSE
)

par(mfrow = c(1, 1))

# ==============================================================================
# 6. Recovering G
# ==============================================================================
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

G <- matrix(0, n_sims, length(y_grid))

cat("\nReconstructing G...\n")

for(i in 1:n_sims) {
  th <- mat_all[i, "theta"]
  si <- mat_all[i, "sigma"]
  mu <- mat_all[i, "mu"]
  
  # Mu fixed at 0
  z <- (y_grid - mu) / si
  
  F_y <- numeric(length(y_grid))
  valid <- z > 0
  F_y[valid] <- plnorm(z[valid], meanlog = mu, sdlog = si)
  F_y <- pmin(pmax(F_y, 1e-15), 1 - 1e-15)
  
  # METHOD A: Exact
  G[i, ] <- gumbel_diag(F_y, th, n_data)
}

# --- True G ---
z_true <- (y_grid - mu_true) / sigma_true
F_true <- numeric(length(y_grid))
valid_true <- z_true > 0
F_true[valid_true] <- plnorm(z_true[valid_true], meanlog = mu_true, sdlog = sigma_true)

G_true <- gumbel_diag(F_true, theta_true, n_data)

# --- Plotting ---
G_mean <- colMeans(G)
G_lower <- apply(G, 2, quantile, 0.05)
G_upper <- apply(G, 2, quantile, 0.95)

par(mfrow = c(1, 1))
plot(y_grid, G_mean, type = "l", lwd = 2, col = "blue",
     ylim = c(0, 1), 
     main = "Estimated Limit Distribution G",
     xlab = "y", ylab = "P(Mn <= y)")

polygon(c(y_grid, rev(y_grid)), c(G_lower, rev(G_upper)), 
        col = rgb(0, 0, 1, 0.1), border = NA)

lines(y_grid, G_true, col = "red", lwd = 2, lty = 2)
