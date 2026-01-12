#' Haario adaptive covariance scheme
#'
#' @param eps Small jitter added to covariance diagonal
#' @param t0 starting adaptation
#'
#' @return Adaptation object

adapt_haario <- function(
  eps = 1e-6,
  t0 = 0
) {
  list(
    type = "haario",
    #' Update proposal covariance using empirical covariance
    #'
    #' @param state Proposal state (contains Sigma, mean, t)
    #' @param param Current parameter vector x_t
    #' @param accept Logical, whether proposal was accepted
    #' @param iter Current iteration number
    update = function(state, param, accept, iter) {
      # note: state values refer to iteration t-1
      # init + first iteration
      state$t <- state$t + 1
      old_mean <- state$mean
      state$mean <- state$mean + (param - state$mean) / state$t # mean update

      # Welford's algorithm
      state$cov <- (state$t - 2) / (state$t - 1) * state$cov + 1/state$t * (tcrossprod(param - old_mean))

      if (iter > t0) {
        d <- length(param)
        state$Sigma <- (2.38^2 / d) * (state$cov + eps * diag(d))
      }

      state
    }
  )
}
