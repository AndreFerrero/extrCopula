# ==============================================================
# t Margin
# ==============================================================
margin_t <- list(

  name = "t",

  # --------------------------
  # Quantile function
  # --------------------------
  simulate = function(u, param) {
    # u: vector of uniforms
    # param: list with mu (location), sigma (scale), df (degrees of freedom)
    qt(u, df = param["df"]) * param["sigma"] + param["mu"]
  },

  # --------------------------
  # Log-density
  # --------------------------
  log_density = function(x, param) {
    dt((x - param["mu"]) / param["sigma"], df = param["df"], log = TRUE) - log(param["sigma"])
  },

  # --------------------------
  # CDF
  # --------------------------
  cdf = function(x, param) {
    pt((x - param["mu"]) / param["sigma"], df = param["df"])
  },

  # --------------------------
  # Log-prior for Bayesian model
  # --------------------------
  log_prior = function(param, prior_mu = c(0, 10), prior_sigma = c(0, 1), prior_df = c(2, 30)) {
    if (param["sigma"] <= 0 || param["df"] <= 2) return(-Inf)
    dnorm(param["mu"], prior_mu[1], prior_mu[2], log = TRUE) +
      dlnorm(param["sigma"], prior_sigma[1], prior_sigma[2], log = TRUE) +
      dunif(param["df"], prior_df[1], prior_df[2], log = TRUE)
  }
)
