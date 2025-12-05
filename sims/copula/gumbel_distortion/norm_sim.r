library(parallel)
library(here)
library(tidyverse)
library(copula)
library(evd)

set.seed(123)

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
n <- 2000            # Sequence length (increased slightly for asymptotic stability)
B <- 500             # Maxima per repetition
R <- 100             # Number of repetitions
alpha <- 2           # Frechet Shape (Tail index)

# We compare Independence vs Weak Dep (1.5) vs Strong Dep (2.5)
scenarios <- c("iid", "theta_1.5", "theta_2.5")

# ==============================================================================
# 2. PARALLEL SETUP
# ==============================================================================
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)

# Load libraries on workers
clusterEvalQ(cl, {
  library(copula)
  library(evd)
})

# Helper function to generate 1 max from Dependent Frechet
# (Defined inline to ensure workers have it)
clusterEvalQ(cl, {
  generate_dep_max <- function(cop, alpha) {
    # 1. Draw n uniform variates from Copula
    u <- rCopula(1, cop) 
    # 2. Transform to Frechet: x = (-log u)^(-1/alpha)
    x <- (-log(u))^(-1/alpha)
    return(max(x))
  }
  
  generate_iid_max <- function(n, alpha) {
    # Standard Frechet generation
    x <- rfrechet(n, shape = alpha)
    return(max(x))
  }
})

clusterExport(cl, c("n", "alpha", "B"), envir = environment())

# ==============================================================================
# 3. SIMULATION LOOP
# ==============================================================================
all_results_by_scenario <- list()

for (scen in scenarios) {
  cat("Running Scenario:", scen, "...\n")
  
  # --- A. DETERMINE NORMALIZING CONSTANTS ---
  if (scen == "iid") {
    theta_val <- 1
    rn <- n
  } else {
    theta_val <- as.numeric(sub("theta_", "", scen))
    rn <- n^(1/theta_val) # The Effective Sample Size
  }
  
  # Theoretical Normalizing Constants for Frechet(alpha)
  # Location b_n = 0
  # Scale a_n = (rn)^(1/alpha)
  bn_theory <- 0
  an_theory <- rn^(1/alpha)
  
  cat("   -> Effective N:", round(rn, 1), "\n")
  cat("   -> Normalizing Scale (an):", round(an_theory, 2), "\n")
  
  # Export specific copula to workers if needed
  if (scen != "iid") {
    clusterExport(cl, "theta_val", envir = environment())
    clusterEvalQ(cl, {
      curr_cop <<- gumbelCopula(param = theta_val, dim = n)
    })
  }
  
  # --- B. GENERATION & FITTING ---
  scenario_results <- vector("list", R)
  
  # Export normalization constants to workers
  clusterExport(cl, c("an_theory", "bn_theory"), envir = environment())
  
  for (r in seq_len(R)) {
    
    # 1. Generate B maxima (Vectorized on workers)
    if (scen == "iid") {
      max_raw <- parSapply(cl, 1:B, function(x) generate_iid_max(n, alpha))
    } else {
      max_raw <- parSapply(cl, 1:B, function(x) generate_dep_max(curr_cop, alpha))
    }
    
    # 2. NORMALIZE
    # Z = (M_n - b_n) / a_n
    # If theory is correct, Z should follow Standard Frechet(alpha) 
    # regardless of theta.
    max_norm <- (max_raw - bn_theory) / an_theory
    
    # 3. Fit GEV to NORMALIZED maxima
    # We expect: loc ~ 1, scale ~ 1/alpha, shape ~ 1/alpha (approx)
    fit <- tryCatch(fgev(max_norm, std.err = FALSE)$estimate,
                    error = function(e) c(NA, NA, NA))
    
    scenario_results[[r]] <- data.frame(
      repetition = r,
      location = fit[1],
      scale    = fit[2],
      shape    = fit[3]
    )
  }
  
  all_results_by_scenario[[scen]] <- bind_rows(scenario_results)
}

stopCluster(cl)

# ==============================================================================
# 4. VISUALIZATION
# ==============================================================================
all_results_long <- bind_rows(lapply(names(all_results_by_scenario), function(scen) {
  df <- all_results_by_scenario[[scen]]
  df$scenario <- scen
  df
}))

all_long <- all_results_long %>%
  pivot_longer(cols = c(location, scale, shape), names_to = "parameter", values_to = "estimate")

# Filter outliers for clean plotting
clean_long <- all_long %>%
  group_by(parameter) %>%
  dplyr::filter(estimate > quantile(estimate, 0.01, na.rm=T) & 
           estimate < quantile(estimate, 0.99, na.rm=T))

# Plot
my_cols <- c("iid" = "#66c2a5", "theta_1.5" = "#fc8d62", "theta_2.5" = "#8da0cb")

p <- ggplot(clean_long, aes(x = scenario, y = estimate, fill = scenario)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~parameter, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = my_cols) +
  labs(title = "GEV Estimates on Normalised Maxima",
       y = "Estimated Parameter") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)   # rotate scenario labels
  )

print(p)

plots_dir <- here("sims", "copula", "gumbel_distortion", "plots")
ggsave(p, file = here(plots_dir, "normalised.pdf"))

# Numeric Summary
summary_stats <- all_results_long %>%
  group_by(scenario) %>%
  summarize(across(c(location, scale, shape), mean, na.rm = TRUE))

print(summary_stats)