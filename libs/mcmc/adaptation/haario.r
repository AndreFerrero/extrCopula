#' Haario adaptive covariance scheme
#'
#' @param eps Small jitter added to covariance diagonal
#' @param t0 starting adaptation
#' 
#' @param check_every check covariance update
#' 
#' @param tol tolderance for covariance check
#' 
#' @param stable_required iterations with covariance under tolerance level
#' 
#' @param t_min minimum adaptation time
#'
#' @return Adaptation object

adapt_haario <- function(
  eps = 1e-6,
  t0 = 0,
  check_every = 100,
  tol = 0.02,
  stable_required = 5,
  t_min = 1000
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

      ready <- FALSE

      if (iter > t0) {
        d <- length(param)
        Sigma_new <- (2.38^2 / d) * (state$cov + eps * diag(d))

        if (iter > t_min && iter %% check_every == 0) {
          rel_change <- norm(Sigma_new - state$Sigma, "F") /
                        norm(state$Sigma, "F")

          state$stable_count <- if (rel_change < tol)
            state$stable_count + 1 else 0

          if (state$stable_count >= stable_required)
            ready <- TRUE
        }

        state$Sigma <- Sigma_new
      }

      list(state = state, ready = ready)
    }
  )
}

