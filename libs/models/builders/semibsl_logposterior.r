# libs/models/builders/semibsl_logposterior.R

source("libs/packages.R")

#' Build a SemiBSL log-posterior function
#'
#' @param simulator Function: simulates n_obs samples given param vector
#' @param sum_stats Function: computes summary statistics from simulated data
#' @param n_sim Number of simulations per likelihood evaluation
#' @param log_prior Function: returns log prior for param
#' @param inverse_transform function to transform parameters in the original space
#' @param log_jacobian jacobian function of the transformation
#'
#' @return Function(param) -> log-posterior
build_semibsl_logposterior <- function(
  copula, margin, param_map,
  data,
  simulator,
  sum_stats,
  n_sim,
  inverse_transform = NULL,
  log_jacobian = NULL,
  margin_prior = NULL,
  copula_prior = NULL
) {

  # Precompute observed summary statistics
  y_obs <- sum_stats(data)

  function(param_init) {

    # Apply inverse transform if given
    # param is always the parameter in the original space
    # no transformation leaves param = param_init
    param <- if (!is.null(inverse_transform)) inverse_transform(param_init) else param_init

    param_m <- param[param_map$margin]
    param_c <- param[param_map$copula]

    # --- Prior ---
    # do.call is to insert the optional prior parameters to overide the defaults
    logprior_m <- do.call(margin$log_prior, c(list(param_m), margin_prior))
    logprior_c <- do.call(copula$log_prior, c(list(param_c), copula_prior))
    logprior <- logprior_m + logprior_c
    
    if (!is.finite(logprior)) return(-Inf)

    # Simulate n_sim datasets and compute summary statistics
    sim_stats <- matrix(NA, nrow = n_sim, ncol = length(y_obs))
    for (i in seq_len(n_sim)) {
      x_sim <- simulator(param, length(data))
      if (any(is.na(x_sim))) return(-Inf)
      sim_stats[i, ] <- sum_stats(x_sim)
    }

    # Marginal KDE log-likelihood
    loglik_marg <- 0
    eta_obs <- numeric(length(y_obs))
    eta_sim <- matrix(NA, nrow = n_sim, ncol = length(y_obs))
    for (j in seq_along(y_obs)) {
      sd_sim <- sd(sim_stats[, j])
      iqr_sim <- IQR(sim_stats[, j])
      h <- 0.9 * min(sd_sim, iqr_sim/1.34) * n_sim^(-0.2)
      if (h == 0) h <- 1e-6

      z_scores <- (y_obs[j] - sim_stats[, j]) / h
      log_vals <- dnorm(z_scores, log = TRUE)
      max_val <- max(log_vals)
      log_pdf <- max_val + log(sum(exp(log_vals - max_val))) - log(n_sim) - log(h)
      loglik_marg <- loglik_marg + log_pdf

      # Gaussian copula transformation
      u_obs <- mean(pnorm(z_scores))
      u_obs <- min(max(u_obs, 1e-6), 1 - 1e-6)
      eta_obs[j] <- qnorm(u_obs)
      u_sim <- rank(sim_stats[, j]) / (n_sim + 1)
      eta_sim[, j] <- qnorm(pmin(pmax(u_sim, 1e-6), 1 - 1e-6))
    }

    # Copula likelihood (Gaussian copula with empirical correlation)
    R_hat <- cor(eta_sim)
    if (det(R_hat) <= 1e-9) R_hat <- R_hat + diag(1e-6, ncol(R_hat))
    R_inv <- tryCatch(solve(R_hat), error = function(e) MASS::ginv(R_hat))
    log_det_R <- as.numeric(determinant(R_hat, logarithm = TRUE)$modulus)
    quad <- t(eta_obs) %*% (R_inv - diag(1, ncol(R_inv))) %*% eta_obs
    loglik_copula <- -0.5 * log_det_R - 0.5 * quad

    # Total semiBSL log-likelihood
    total_loglik <- loglik_marg + loglik_copula

    # Jacobian adjustment if in transformed space
    logjac <- if (!is.null(log_jacobian)) log_jacobian(param_init) else 0

    logprior + total_loglik + logjac
  }
}
