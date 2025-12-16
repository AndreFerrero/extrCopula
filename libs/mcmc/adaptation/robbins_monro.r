#' Robbins–Monro adaptation of proposal scale
#'
#' Adapts a scalar scale factor to target a desired acceptance rate.
#'
#' @param target_accept Target acceptance probability (default 0.234)
#' @param gamma0 Initial step size
#' @param kappa Decay exponent (must be > 0.5)
#'
#' @return Adaptation object
adapt_robbins_monro <- function(
  target_accept = 0.234,
  gamma0 = 0.05,
  kappa = 0.6
) {

  list(

    update = function(state, param, accept, iter) {

      if (iter == 1) {
        if (is.null(state$scale)) state$scale <- 1
        if (is.null(state$Sigma0)) state$Sigma0 <- state$Sigma
        return(state)
      }

      # Robbins–Monro step size
      gamma_t <- gamma0 / iter^kappa

      # Log-scale update
      log_scale <- log(state$scale) +
        gamma_t * (as.numeric(accept) - target_accept)

      state$scale <- exp(log_scale)

      # Update proposal covariance
      state$Sigma <- state$scale^2 * state$Sigma0

      state
    }
  )
}
