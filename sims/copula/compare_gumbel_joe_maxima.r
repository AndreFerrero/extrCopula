library(ggplot2)
library(dplyr)

source("libs/models/copulas/gumbel.R")     
source("libs/models/copulas/joe.r")    
source("libs/models/margins/lognormal.R")
source("libs/models/margins/t.R")
source("libs/models/margins/frechet.R")
source("libs/models/builders/simulator.R") 

# ------------------------------
# 0. Parameter map (shared)
# ------------------------------
param_map <- list(
  margin = c("mu", "sigma", "df", "scale", "shape")[c(TRUE, TRUE, TRUE, TRUE, TRUE)],
  copula = "theta"
)

# ------------------------------
# 1. Define simulation function
# ------------------------------
simulate_maxima <- function(copula_obj, margin_obj, param, n_obs, n_reps) {
  maxima <- numeric(n_reps)
  simulator <- build_simulator(copula_obj, margin_obj, param_map)
  
  for (i in 1:n_reps) {
    X <- simulator(param, n_obs)
    if (all(is.na(X))) stop("Simulator returned NA")
    maxima[i] <- max(X, na.rm = TRUE)
  }
  
  maxima
}

# ------------------------------
# 2. Copula and margin combos
# ------------------------------
copula_list <- list(
  gumbel = copula_gumbel,
  joe    = copula_joe
)

margin_list <- list(
  lognormal = margin_lognormal,
  frechet   = margin_frechet,
  t         = margin_t
)

theta_vec <- c(1, 2, 3)      # copula dependence
n_obs <- 1000
n_reps <- 200

# Define parameter sets per margin
param_list <- list(
  lognormal = c(mu = 0, sigma = 1, theta = 2),
  t         = c(mu = 0, sigma = 1, df = 3, theta = 2),
  frechet   = c(scale = 1, shape = 2, theta = 2)
)

results <- list()

# ------------------------------
# 3. Simulation loop
# ------------------------------
for (margin_name in names(margin_list)) {
  margin_obj <- margin_list[[margin_name]]
  base_param <- param_list[[margin_name]]
  
  for (cop_name in names(copula_list)) {
    copula_obj <- copula_list[[cop_name]]
    
    for (th in theta_vec) {
      param <- base_param
      param["theta"] <- th
      cat("Simulating for", cop_name, "theta =", th, "with margin", margin_name, "...\n")
      maxima <- simulate_maxima(copula_obj, margin_obj, param, n_obs, n_reps)
      
      results[[paste(cop_name, margin_name, th, sep = "_")]] <- data.frame(
        copula = cop_name,
        margin = margin_name,
        theta  = th,
        maxima = maxima
      )
    }
  }
}

maxima_df <- bind_rows(results)

# ------------------------------
# 4. Plot maxima densities
# ------------------------------
maxima_df$theta <- factor(maxima_df$theta, levels = theta_vec)

ggplot(maxima_df, aes(x = maxima, color = theta, linetype = copula)) +
  geom_density(size = 1.2) +
  scale_x_log10() +  # logarithmic x-axis
  facet_wrap(~margin, scales = "free") +
  labs(
    title = "Effect of Copula Dependence and Margin on Maxima (Log Scale)",
    x = "Maxima (log scale)", y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# ------------------------------
# 5. Summarize key quantiles
# ------------------------------
summary_df <- maxima_df %>%
  group_by(copula, margin, theta) %>%
  summarize(
    median_max = median(maxima),
    q90 = quantile(maxima, 0.9),
    q95 = quantile(maxima, 0.95),
    q99 = quantile(maxima, 0.99),
    .groups = "drop"
  )

print(summary_df)
