# ------------------------------------------------------------
# Benchmark simulation cost vs number of observations
# ------------------------------------------------------------

source("libs/packages.R")
source("libs/models/copulas/gumbel.R")
source("libs/models/margins/lognormal.R")
source("libs/models/builders/simulator.R")

param_map <- list(margin = c("mu", "sigma"), copula = "theta")
simulator <- build_simulator(copula_gumbel, margin_lognormal, param_map)

set.seed(123)
true_param <- c(mu = 0, sigma = 1, theta = 2)

# Observation sizes to test
n_grid <- c( 500, 1000, 2000, 5000, 10000, 15000, 20000)

# Number of repetitions per n (to smooth timing noise)
n_rep <- 100

timings <- data.frame(
  n_obs = n_grid,
  mean_time = NA,
  sd_time = NA
)

for (k in seq_along(n_grid)) {

  n_obs <- n_grid[k]
  times <- numeric(n_rep)

  for (r in seq_len(n_rep)) {
    t0 <- proc.time()[3]
    x <- simulator(true_param, n_obs)
    times[r] <- proc.time()[3] - t0
  }

  timings$mean_time[k] <- mean(times)
  timings$sd_time[k]   <- sd(times)

  cat(sprintf("n = %6d | mean %.4f s\n", n_obs, timings$mean_time[k]))
}

timings
log(timings)
plot(
  timings$n_obs,
  timings$mean_time,
  type = "b",
  pch = 19,
  xlab = "Number of observations (n)",
  ylab = "Simulation time (seconds)",
  main = "Simulation cost vs number of observations"
)

plot(
  log(timings$n_obs),
  log(timings$mean_time),
  type = "b",
  pch = 19,
  xlab = "log(n)",
  ylab = "log(time)",
  main = "Log–log scaling of simulation cost"
)
abline(lm(log(mean_time) ~ log(n_obs), data = timings), col = "red")

lm(log(mean_time) ~ log(n_obs), data = timings)
