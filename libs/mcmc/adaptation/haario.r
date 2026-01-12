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
    #' @param param Current parameter vector
    #' @param accept Logical, whether proposal was accepted
    #' @param iter Current iteration number
    update = function(state, param, accept, iter) {
      state$t <- state$t + 1
      delta <- param - state$mean
      state$mean <- state$mean + delta / state$t
      state$cov <- state$cov + (tcrossprod(delta) - state$cov) / state$t

      if (iter > t0) {
        d <- length(param)
        state$Sigma <- (2.38^2 / d) * (state$cov + eps * diag(d))
      }

      list(state = state)
    }
  )
}
