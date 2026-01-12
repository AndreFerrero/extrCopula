#' Gaussian random-walk proposal
#'
#' @param Sigma0 Initial proposal covariance matrix (p x p)
#'
#' @return A proposal object with:
#'   - init_state(param): initializes internal proposal state
#'   - propose(param, state): proposes new parameter vector
proposal_gaussian_rw <- function(Sigma0) {

  list(
    #' Initialize proposal state
    #' @param param Current parameter vector
    #' @return List storing proposal state
    init_state = function(param) {
      list(
        Sigma  = Sigma0,              # current covariance used for proposals
        Sigma0 = Sigma0,              # reference covariance for scaling
        cov    = Sigma0
,              # empirical covariance
        mean   = param,              # running mean of chain
        t      = 1,                  # iteration counter
        scale  = 1                   # Robbins–Monro scale
      )
    },

    #' Propose a new parameter value
    #' @param param Current parameter vector
    #' @param state Current proposal state
    #' @return List with:
    #'   - param: proposed parameter vector
    #'   - log_q_ratio: log symmetric proposal ratio (always 0)
    propose = function(param, state) {
      eps <- MASS::mvrnorm(
        n = 1,
        mu = rep(0, length(param)),
        Sigma = state$Sigma
      )
      list(param = param + eps, log_q_ratio = 0)
    }
  )
}
