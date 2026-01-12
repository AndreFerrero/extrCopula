#' Hybrid adaptive scheme: Robbins–Monro + Haario covariance
#'
#' @param target_accept Target acceptance probability for Robbins–Monro
#' @param gamma0 Initial Robbins–Monro step size
#' @param kappa Robbins–Monro decay exponent (> 0.5)
#' @param eps Small jitter added to covariance diagonal
#' @param t0 Start adapting empirical covariance (Haario)
#'
#' @return Adaptation object with update(state, param, accept, iter)

adapt_hybrid <- function(
  target_accept = 0.234,
  gamma0 = 1,
  kappa = 0.6, # 0.6 is often more stable than 1
  eps = 1e-6,
  t0 = 0
) {
  list(
    type = "hybrid",
    update = function(state, param, accept, iter) {
      d <- length(param)

      # --- 1. Robbins–Monro: adapt scale (runs from the start) ---
      # Targetting the scale of the jumps
      gamma_t <- gamma0 / iter^kappa
      log_scale <- log(state$scale) + gamma_t * (as.numeric(accept) - target_accept)
      state$scale <- exp(log_scale)

      # note: state values refer to iteration t-1
      # init + first iteration
      state$t <- state$t + 1
      old_mean <- state$mean
      state$mean <- state$mean + (param - state$mean) / state$t # mean update

      # Welford's algorithm
      state$cov <- (state$t - 2) / (state$t - 1) * state$cov + 1/state$t * (tcrossprod(param - old_mean))

      # --- 2. Haario: adapt shape (starts after t0) ---
      if (iter > t0) {
        # Base matrix is now the empirical covariance
        base_matrix <- state$cov
      } else {
        # Before t0, use the user-provided Sigma0 as the 'shape'
        base_matrix <- state$Sigma0
      }

      # --- 3. Update proposal covariance ---
      # Paper convention: Sigma = lambda * Cov
      # Using scale^2 here matches your logic
      state$Sigma <- state$scale^2 * (base_matrix + eps * diag(d))

      state
    }
  )
}
