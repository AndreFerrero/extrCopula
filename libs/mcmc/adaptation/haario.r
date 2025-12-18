#' Haario adaptive covariance scheme
#'
#' @param eps Small jitter added to covariance diagonal
#' @param t0 stating adaptation
#' 
#' @return Adaptation object
adapt_haario <- function(eps = 1e-6, t0 = 1000) {

  list(

    #' Update proposal covariance using empirical covariance
    #'
    #' @param state Proposal state (contains Sigma, mean, t)
    #' @param param Current parameter vector
    #' @param accept Logical, whether proposal was accepted
    #' @param iter Current iteration number
    update = function(state, param, accept, iter) {

      if (iter <= t0) return(state)

      if (iter == t0 + 1) {
        state$mean <- param
        state$cov  <- diag(1e-3, length(param))
        state$t <- 1
        return(state)
      }

      # online updating of mean and covariance to avoid storing chain
      state$t <- state$t + 1

      delta <- param - state$mean
      state$mean <- state$mean + delta / state$t

      state$cov <- state$cov +
        (tcrossprod(delta) - state$cov) / state$t

      d <- length(param)
      state$Sigma <- (2.38^2 / d) * (state$cov + eps * diag(d))

      state
    }
  )
}
