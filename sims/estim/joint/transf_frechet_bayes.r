library(copula)
library(evd)
library(coda)
library(mvtnorm)
library(future)
library(future.apply)
library(psych)

set.seed(123)

# -----------------------------------------------------------
# 1. Simulate Data (Fréchet Margins)
# -----------------------------------------------------------
n <- 150
theta_true <- 3
alpha_true <- 3
mu_true <- 10
sigma_true <- 5

cop <- gumbelCopula(param = theta_true, dim = n)
U <- as.numeric(rCopula(1, cop))
X <- mu_true + sigma_true * (-log(U))^(-1 / alpha_true)

# -----------------------------------------------------------
# 2. Fully Transformed Posterior
# The problem of this approach is that transforming mu makes the likelihood unstable when mu -> min(x)
# -----------------------------------------------------------

log_posterior_full_transform <- function(param_vec, data, data_min) {
    # 1. UNPACK TRANSFORMED PARAMETERS
    lambda <- param_vec[1] # for theta
    phi <- param_vec[2] # for mu
    psi <- param_vec[3] # for sigma (log-sigma)
    tau <- param_vec[4] # for alpha (log-alpha)

    # 2. INVERSE TRANSFORMS (Map to Natural Scale)
    theta <- exp(lambda) + 1
    mu <- data_min - exp(phi)
    sigma <- exp(psi)
    alpha <- exp(tau)

    # 3. PRIORS (Defined on Natural Scale)
    if (theta > 20) {
        return(-Inf)
    }
    lp_theta <- dgamma(theta - 1, shape = 2, rate = 1, log = TRUE)

    # Sigma: Boundary Avoiding
    lp_sigma <- dlnorm(sigma, meanlog = log(sd(data) / 2), sdlog = 0.5, log = TRUE)

    # Alpha: Gamma
    lp_alpha <- dgamma(alpha, shape = 3, rate = 1, log = TRUE)

    # Mu: Normal
    lp_mu <- dnorm(mu, mean = mean(data) - sd(data), sd = 10, log = TRUE)

    # 4. LIKELIHOODS
    # Marginal
    log_dens_vals <- dfrechet(data, mu, sigma, alpha, log = TRUE)
    if (any(log_dens_vals == -Inf)) {
        return(-Inf)
    }

    ll_margin <- sum(log_dens_vals)

    # Copula
    u_hat <- pfrechet(data, mu, sigma, alpha)
    u_hat <- pmin(pmax(u_hat, 1e-8), 1 - 1e-8)

    ld_cop <- dCopula(u_hat,
        copula = gumbelCopula(theta, dim = length(data)),
        log = TRUE
    )

    if (is.na(ld_cop) || !is.finite(ld_cop)) {
        return(-Inf)
    }

    # 5. JACOBIAN ADJUSTMENT
    log_jacobian <- lambda + phi + psi + tau

    return(lp_theta + lp_mu + lp_sigma + lp_alpha +
        ll_margin + ld_cop +
        log_jacobian)
}

# -----------------------------------------------------------
# 3. Worker (Unconstrained Space)
# -----------------------------------------------------------
run_worker <- function(data, n_iter, init_transformed, chain_id, progress_dir) {
    library(copula)
    library(mvtnorm)
    library(coda)
    library(evd)

    p_names <- c("lambda", "phi", "psi", "tau")
    n_par <- 4
    chain <- matrix(NA, nrow = n_iter, ncol = n_par)
    colnames(chain) <- p_names

    data_min <- min(data)
    curr_par <- init_transformed

    curr_lp <- log_posterior_full_transform(curr_par, data, data_min)
    if (!is.finite(curr_lp)) stop("Invalid Start")

    # Adaptive Settings
    cov_mat <- diag(rep(0.01, n_par))
    sd_scale <- (2.38^2) / n_par
    eps <- 1e-6 * diag(n_par)
    adapt_start <- 500
    total_accepts <- 0
    status_file <- file.path(progress_dir, paste0("status_", chain_id, ".txt"))

    for (i in 1:n_iter) {
        # Propose in Transformed Space
        prop_val <- rmvnorm(1, mean = curr_par, sigma = cov_mat)[1, ]
        prop_lp <- log_posterior_full_transform(prop_val, data, data_min)
        ratio <- prop_lp - curr_lp

        if (is.finite(ratio) && log(runif(1)) < ratio) {
            curr_par <- prop_val
            curr_lp <- prop_lp
            total_accepts <- total_accepts + 1
        }
        chain[i, ] <- curr_par

        # Adapt
        if (i > adapt_start && i %% 50 == 0) {
            emp_cov <- cov(chain[1:i, ])
            cov_mat <- (sd_scale * emp_cov) + eps
        }

        # Write Status
        if (i %% 50 == 0 || i == n_iter) {
            try(
                {
                    writeLines(sprintf("%d|%.1f", as.integer(i / n_iter * 100), total_accepts / i * 100), status_file)
                },
                silent = TRUE
            )
        }
    }
    return(mcmc(chain))
}

# -----------------------------------------------------------
# 4. Execution
# -----------------------------------------------------------
n_chains <- 3
plan(multisession, workers = n_chains)
progress_dir <- tempdir()
file.remove(list.files(progress_dir, pattern = "status_", full.names = TRUE))

n_iter <- 15000

# Initialize in TRANSFORMED SPACE
inits_list <- list()
d_min <- min(X)

for (c in 1:n_chains) {
    th_g <- runif(1, 1.5, 6)
    mu_g <- d_min - runif(1, 0.5, 4)
    si_g <- runif(1, 2, 7)
    al_g <- runif(1, 2, 7)

    inits_list[[c]] <- c(
        lambda = log(th_g - 1),
        phi    = log(d_min - mu_g),
        psi    = log(si_g),
        tau    = log(al_g)
    )
}

cat("MCMC Started (Fully Reparameterized)...\n")
futures_list <- list()
for (c in 1:n_chains) {
    futures_list[[c]] <- future(
        {
            run_worker(X, n_iter, inits_list[[c]], c, progress_dir)
        },
        seed = TRUE
    )
}

# --- ROBUST PROGRESS BAR ---
last_print_time <- Sys.time()
print_interval <- 2

while (!all(resolved(futures_list))) {
    if (difftime(Sys.time(), last_print_time, units = "secs") > print_interval) {
        status_texts <- character(n_chains)
        for (c in 1:n_chains) {
            fpath <- file.path(progress_dir, paste0("status_", c, ".txt"))
            if (file.exists(fpath)) {
                info <- suppressWarnings(try(readLines(fpath, n = 1), silent = TRUE))
                if (inherits(info, "try-error") || length(info) == 0) {
                    status_texts[c] <- "[Init]"
                } else {
                    parts <- strsplit(info, "\\|")[[1]]
                    if (length(parts) >= 2) {
                        status_texts[c] <- sprintf("C%d: %3s%% (Acc: %4s%%)", c, parts[1], parts[2])
                    } else {
                        status_texts[c] <- "[Err]"
                    }
                }
            } else {
                status_texts[c] <- "[Pending]"
            }
        }
        cat(format(Sys.time(), "%H:%M:%S"), "|", paste(status_texts, collapse = " | "), "\n")
        flush.console()
        last_print_time <- Sys.time()
    }
    Sys.sleep(0.5)
}
cat("Done.\n")

# -----------------------------------------------------------
# 5. Post-Processing
# -----------------------------------------------------------
res_dir <- here("sims", "estim", "joint", "res")

save(futures_list, file = here(res_dir, "transf_frechet_bayes_chains.Rdata"))
load(here(res_dir, "trans_frechet_bayes_chains.Rdata"))

chains_raw <- value(futures_list)

# Explicitly pass d_min to the environment for inversion
chains_inv <- lapply(chains_raw, function(ch) {
    theta <- exp(ch[, 1]) + 1
    mu <- d_min - exp(ch[, 2]) # d_min is from global scope
    sigma <- exp(ch[, 3])
    alpha <- exp(ch[, 4])
    cbind(theta, mu, sigma, alpha)
})

for (i in 1:length(chains_inv)) colnames(chains_inv[[i]]) <- c("theta", "mu", "sigma", "alpha")

burn_in <- n_iter / 2
mcmc_clean <- window(mcmc.list(lapply(chains_inv, mcmc)), start = burn_in + 1, thin = 5)

cat("\n--- DIAGNOSTICS ---\n")
print(gelman.diag(mcmc_clean))
print(effectiveSize(mcmc_clean))

# Visuals
par(mfrow = c(2, 2))
plot(mcmc_clean[, "theta"], main = "Theta", density = FALSE, auto.layout = FALSE)
plot(mcmc_clean[, "mu"], main = "Mu", density = FALSE, auto.layout = FALSE)
plot(mcmc_clean[, "sigma"], main = "Sigma", density = FALSE, auto.layout = FALSE)
plot(mcmc_clean[, "alpha"], main = "Alpha", density = FALSE, auto.layout = FALSE)
par(mfrow = c(1, 1))

# Pairs Plot
mat_all <- as.matrix(mcmc_clean)
pairs.panels(mat_all, method = "pearson", hist.col = "skyblue", density = TRUE, ellipses = TRUE)

cat("\nTrue Values: Theta=3, Mu=10, Sigma=5, Alpha=3\n")
cat("Posterior Means:\n")
print(colMeans(mat_all))
