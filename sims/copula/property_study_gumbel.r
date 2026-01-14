# ------------------------------
# 1. Load required packages
# ------------------------------
source("libs/packages.R")

# ------------------------------
# 2. Source your libs
# ------------------------------
source("libs/models/copulas/gumbel.R")
source("libs/models/copulas/clayton.r")
source("libs/models/margins/lognormal.R")
source("libs/models/builders/simulator.R")

# ------------------------------
# 3. Define parameter mapping
# ------------------------------
param_map <- list(
  margin = c("mu", "sigma"),
  copula = "theta"
)

# ------------------------------
# 4. Build the simulator
# ------------------------------
simulator <- build_simulator(
  copula = copula_gumbel,
  margin = margin_lognormal,
  param_map = param_map
)

# ============================================================
# STUDY OF DEPENDENCE IN EXCHANGEABLE GUMBEL–MARGIN MODEL
# ============================================================

set.seed(123)

# ------------------------------------------------------------
# Study setup
# ------------------------------------------------------------
n <- 1000
n_rep <- 1000
theta_vec <- c(1, 1.5, 2, 3, 5)
param_base <- c(mu = 0, sigma = 1)
u_prob <- 0.95

# Storage
results <- list()

# ------------------------------------------------------------
# Functions for each property
# ------------------------------------------------------------
order_stat_study <- function(theta) {
  max_vals <- numeric(n_rep)
  spacing_vals <- numeric(n_rep)
  for (r in 1:n_rep) {
    X <- simulator(c(param_base, theta = theta), n)
    X_ord <- sort(X)
    max_vals[r] <- X_ord[n]
    spacing_vals[r] <- X_ord[n] - X_ord[n-1]
  }
  list(
    max_mean = mean(max_vals),
    max_sd = sd(max_vals),
    spacing_mean = mean(spacing_vals),
    spacing_sd = sd(spacing_vals)
  )
}

exceedance_count_study <- function(theta) {
  counts <- numeric(n_rep)
  u <- qlnorm(u_prob, meanlog = param_base["mu"], sdlog = param_base["sigma"])
  for (r in 1:n_rep) {
    X <- simulator(c(param_base, theta = theta), n)
    counts[r] <- sum(X > u)
  }
  list(
    mean_count = mean(counts),
    var_count = var(counts),
    dispersion = var(counts) / mean(counts),
    sd_count = sd(counts)
  )
}

conditional_exceedance_study <- function(theta) {
  u <- qlnorm(u_prob, meanlog = param_base["mu"], sdlog = param_base["sigma"])
  cond_vals <- numeric(n_rep)
  unconditional_vals <- numeric(n_rep)
  for (r in 1:n_rep) {
    X <- simulator(c(param_base, theta = theta), n)
    exceed <- which(X > u)
    unconditional_vals[r] <- mean(X > u)
    if (length(exceed) >= 2) {
      cond_vals[r] <- (length(exceed)-1)/(n-1)
    } else {
      cond_vals[r] <- NA
    }
  }
  list(
    unconditional = mean(unconditional_vals),
    conditional = mean(cond_vals, na.rm = TRUE),
    conditional_sd = sd(cond_vals, na.rm = TRUE)
  )
}

edf_variance_study <- function(theta, x0 = 1) {
  Fhat_vals <- numeric(n_rep)
  for (r in 1:n_rep) {
    X <- simulator(c(param_base, theta = theta), n)
    Fhat_vals[r] <- mean(X <= x0)
  }
  list(
    mean_Fhat = mean(Fhat_vals),
    var_Fhat = var(Fhat_vals),
    sd_Fhat = sd(Fhat_vals)
  )
}

summary_variability_study <- function(theta) {
  means <- numeric(n_rep)
  q95 <- numeric(n_rep)
  for (r in 1:n_rep) {
    X <- simulator(c(param_base, theta = theta), n)
    means[r] <- mean(X)
    q95[r] <- quantile(X, 0.95)
  }
  list(
    mean_of_means = mean(means),
    var_of_means = var(means),
    sd_of_means = sd(means),
    mean_q95 = mean(q95),
    var_q95 = var(q95),
    sd_q95 = sd(q95)
  )
}

# ------------------------------------------------------------
# Run full study
# ------------------------------------------------------------
for (theta in theta_vec) {
  cat("Running theta =", theta, "\n")
  results[[as.character(theta)]] <- list(
    order_stats = order_stat_study(theta),
    exceedances = exceedance_count_study(theta),
    conditional = conditional_exceedance_study(theta),
    edf_var = edf_variance_study(theta),
    summaries = summary_variability_study(theta)
  )
}

# ------------------------------------------------------------
# Prepare comparison data frame with SDs for plotting
# ------------------------------------------------------------
comparison <- data.frame(
  theta = theta_vec,
  max_mean = sapply(results, function(x) x$order_stats$max_mean),
  max_sd = sapply(results, function(x) x$order_stats$max_sd),
  spacing_mean = sapply(results, function(x) x$order_stats$spacing_mean),
  spacing_sd = sapply(results, function(x) x$order_stats$spacing_sd),
  dispersion_mean = sapply(results, function(x) x$exceedances$dispersion),
  dispersion_sd = sapply(results, function(x) x$exceedances$sd_count),
  cond_mean = sapply(results, function(x) x$conditional$conditional),
  cond_sd = sapply(results, function(x) x$conditional$conditional_sd),
  var_mean = sapply(results, function(x) x$summaries$var_of_means),
  var_sd = sapply(results, function(x) x$summaries$sd_of_means)
)

print(comparison)

# ------------------------------------------------------------
# Plot statistics with uncertainty bands
# ------------------------------------------------------------
library(ggplot2)
library(tidyr)
# Select only the main statistics for plotting
plot_stats <- comparison %>%
  select(theta, max_mean, spacing_mean, dispersion_mean, cond_mean, var_mean)

# Pivot for ggplot
plot_df <- plot_stats %>%
  pivot_longer(cols = -theta, names_to = "stat", values_to = "value") %>%
  mutate(
    lower = case_when(
      stat == "max_mean"     ~ value - comparison$max_sd,
      stat == "spacing_mean" ~ value - comparison$spacing_sd,
      stat == "dispersion_mean" ~ pmax(value - comparison$dispersion_sd, 0),
      stat == "cond_mean"       ~ pmax(value - comparison$cond_sd, 0),
      stat == "var_mean"        ~ pmax(value - comparison$var_sd, 0)
    ),
    upper = case_when(
      stat == "max_mean"     ~ value + comparison$max_sd,
      stat == "spacing_mean" ~ value + comparison$spacing_sd,
      stat == "dispersion_mean" ~ value + comparison$dispersion_sd,
      stat == "cond_mean"       ~ pmin(value + comparison$cond_sd, 1),
      stat == "var_mean"        ~ value + comparison$var_sd
    )
  )

# Rename stats
stat_labels <- c(
  max_mean = "Mean Maximum",
  spacing_mean = "Mean Upper Spacing",
  dispersion_mean = "Dispersion of Exceedances",
  cond_mean = "Conditional Exceedance Probability",
  var_mean = "Variance of Sample Mean"
)
plot_df$stat_label <- stat_labels[plot_df$stat]

# Plot
ggplot(plot_df, aes(x = theta, y = value)) +
  geom_line(aes(color = stat_label), size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = stat_label), alpha = 0.2) +
  geom_point(aes(color = stat_label), size = 2) +
  facet_wrap(~stat_label, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 14) +
  labs(
    x = expression(theta),
    y = "Statistic Value",
    title = "Effect of Theta on Key Sample Statistics with Uncertainty Bands",
    color = "Statistic",
    fill = "Statistic"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12)
  )

