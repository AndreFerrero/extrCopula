#' Run a single MCMC chain
#'
#' @param log_target Function returning log post(param)
#' @param init Initial parameter vector
#' @param n_iter Number of MCMC iterations
#' @param proposal Proposal object
#' @param adapt Adaptation object
#'
#' @return List with samples matrix and acceptance rate

source("libs/packages.R")

run_chain <- function(
  log_target,
  init,
  n_iter,
  proposal,
  adapt = adapt_none()
) {

  p <- length(init)

  param <- init
  logpost  <- log_target(param)

  samples <- matrix(NA, n_iter, p)
  accept  <- logical(n_iter)

  # Initialize proposal state
  prop_state <- proposal$init_state(param)

  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)

  for (i in seq_len(n_iter)) {

    # MH step, returns proposal with acceptance and logposterior info
    step <- mh_step(
      param = param,
      logpost = logpost,
      log_target = log_target,
      proposal = proposal,
      prop_state = prop_state
    )

    param <- step$param
    logpost  <- step$logpost
    accept[i] <- step$accept

    # Adapt proposal
    prop_state <- adapt$update(
      prop_state, param, accept[i], i
    )

    samples[i, ] <- param

    # Update progress bar
    setTxtProgressBar(pb, i)
  }

  close(pb)
  
  colnames(samples) <- names(init)

  list(
    samples = samples,
    accept_rate = mean(accept),
    conv_state = prop_state
  )
}
