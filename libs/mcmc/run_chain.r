run_chain <- function(
  log_target,
  init,
  n_iter,
  proposal,
  burn_in,
  engine_step = mh_step,
  adapt = adapt_none()
) {
  p <- length(init)

  param <- init
  logpost <- log_target(param)

  samples <- matrix(NA, n_iter, p)
  accept <- logical(n_iter)

  # Initialize proposal state
  prop_state <- proposal$init_state(param)

  # Adaptation flag
  adapting <- adapt$type != "none" && burn_in > 0

  # Progress bar
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)

  for (i in seq_len(n_iter)) {
    # MH step
    step <- engine_step(
      param      = param,
      logpost    = logpost,
      log_target = log_target,
      proposal   = proposal,
      prop_state = prop_state
    )

    param <- step$param
    logpost <- step$logpost
    accept[i] <- step$accept
    samples[i, ] <- param

    # Adaptation during burn_in
    if (adapting && i <= burn_in) {
      res <- adapt$update(prop_state, param, accept[i], i)
      prop_state <- res$state
    }

    # Force adaptation stop at burn-in
    if (adapting && i == burn_in) {
      adapting <- FALSE
      adapt <- adapt_none()
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)
  colnames(samples) <- names(init)

  list(
    samples      = samples,
    burn_in       = burn_in,
    accept_rate  = mean(accept),
    conv_state   = prop_state
  )
}
