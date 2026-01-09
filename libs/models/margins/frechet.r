# ==============================================================
# Frechet Margin
# ==============================================================
margin_frechet <- list(

  name = "frechet",

  # --------------------------
  # Quantile function (inverse CDF)
  # --------------------------
  simulate = function(u, param) {
    # param: scale = s > 0, shape = alpha > 0
    s <- param["scale"]
    alpha <- param["shape"]
    if (s <= 0 || alpha <= 0) stop("Invalid Frechet parameters")
    s * (-log(u))^(-1 / alpha)
  },

  # --------------------------
  # Log-density
  # --------------------------
  log_density = function(x, param) {
    s <- param["scale"]
    alpha <- param["shape"]
    if (any(x <= 0) || s <= 0 || alpha <= 0) return(-Inf)
    log(alpha) - alpha * log(s) - (1 + alpha) * log(x) - (s / x)^alpha
  },

  # --------------------------
  # CDF
  # --------------------------
  cdf = function(x, param) {
    s <- param["scale"]
    alpha <- param["shape"]
    if (any(x <= 0)) return(0)
    exp(-(s / x)^alpha)
  },

  # --------------------------
  # Log-prior for Bayesian model
  # --------------------------
  log_prior = function(param, prior_s = c(0.1, 10), prior_alpha = c(0.1, 5)) {
    s <- param["scale"]
    alpha <- param["shape"]
    if (s <= 0 || alpha <= 0) return(-Inf)
    dgamma(s, shape = 2, rate = 1, log = TRUE) +
      dgamma(alpha, shape = 2, rate = 1, log = TRUE)
  }
)
