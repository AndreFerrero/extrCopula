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
      prop_state <- adapt$update(prop_state, param, accept[i], i)
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

  burn_in_accept_rate <- if (burn_in > 0) {
    mean(accept[1:burn_in])
  } else {
    NA_real_
  }

  after_burn_in_accept_rate <- if (burn_in < n_iter) {
    mean(accept[(burn_in + 1):n_iter])
  } else {
    NA_real_
  }

  list(
    samples = samples,
    burn_in = burn_in,
    accept_rate = mean(accept),
    after_burn_in_accept_rate = after_burn_in_accept_rate,
    burn_in_accept_rate = burn_in_accept_rate,
    conv_state = prop_state
  )
}
