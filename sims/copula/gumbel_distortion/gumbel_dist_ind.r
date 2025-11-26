library(parallel)
library(here)
library(tidyverse)
library(copula)
library(evd)

set.seed(123)

plots_dir <- here("sims", "copula", "gumbel_distortion", "plots")
res_dir   <- here("sims", "copula", "gumbel_distortion", "res")
common_dir <- here("sims", "common")

# ============================================================
# Helper functions
# ============================================================
source(here(common_dir, "handy_funs.r"))   # must provide rCopFrechet etc.

# ============================================================
# Simulation settings (change as needed)
# ============================================================
n <- 20000            # length of each sequence X_1...X_n
B <- 1000             # how many maxima per repetition (M_1 ... M_B)
R <- 100              # how many repetitions of the whole experiment
alpha <- 2
scenarios <- c("iid", "theta_1.5", "theta_2.5")

# ============================================================
# Parallel setup
# ============================================================
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {
  library(copula)
  library(evd)
  # ensure helper functions available on nodes if needed
})
# Export data/funcs commonly needed on workers
clusterExport(cl, c("n", "alpha", "B", "rCopFrechet"), envir = environment())

# ============================================================
# Simulation: for each scenario, repeat R times:
#   - for each repetition: generate B maxima (each max from an independent sample of size n)
#   - fit GEV (fgev) to the B maxima -> get estimates (location, scale, shape)
# ============================================================
# container for results
# structure: list[[scenario]] -> data.frame with columns: repetition, location, scale, shape
all_results_by_scenario <- list()

for (scen in scenarios) {
  cat("Starting scenario:", scen, "...\n")
  scenario_results <- vector("list", R)

  if (scen == "iid") {
    # iid Frechet(alpha) sequences
    for (r in seq_len(R)) {
      # generate B maxima in parallel: for each b, draw n frechet variates, take max
      max_vec <- parSapply(cl, 1:B, function(i, n, alpha) {
        X <- rfrechet(n, shape = alpha)    # from evd
        max(X)
      }, n, alpha)

      # fit GEV to the B maxima
      fit <- tryCatch(fgev(max_vec, std.err = FALSE)$estimate,
                      error = function(e) c(NA, NA, NA))
      scenario_results[[r]] <- data.frame(
        repetition = r,
        location = fit[1],
        scale    = fit[2],
        shape    = fit[3]
      )
      if (r %% 10 == 0) cat("  repetition", r, "done\n")
    }

  } else {
    # Gumbel copula case: extract numeric theta from scen
    tval <- as.numeric(sub("theta_", "", scen))
    # Export tval so worker sees it
    clusterExport(cl, c("tval"), envir = environment())
    # construct a Gumbel copula object on the workers (we want to sample n-dim copula)
    # note: many copula implementations expect dimension to match the vector length
    clusterEvalQ(cl, {
      # Gcop will be used inside worker function; set dimension to n
      Gcop <<- gumbelCopula(param = tval, dim = n)
    })

    for (r in seq_len(R)) {
      # For each b: draw n-dim sample from copula and transform marginal to Frechet via rCopFrechet
      max_vec <- parSapply(cl, 1:B, function(i, alpha) {
        # rCopFrechet is assumed to take (alpha, Gcop) and produce a length-n Frechet vector
        X <- rCopFrechet(alpha, Gcop)
        max(X)
      }, alpha)

      fit <- tryCatch(fgev(max_vec, std.err = FALSE)$estimate,
                      error = function(e) c(NA, NA, NA))
      scenario_results[[r]] <- data.frame(
        repetition = r,
        location = fit[1],
        scale    = fit[2],
        shape    = fit[3]
      )
      if (r %% 10 == 0) cat("  repetition", r, "done\n")
    }

    # clean Gcop on workers if you want:
    clusterEvalQ(cl, { if (exists("Gcop")) rm(Gcop) })
  }

  all_results_by_scenario[[scen]] <- bind_rows(scenario_results)
}

stopCluster(cl)

# ============================================================
# Combine into one tidy data frame for analysis
# ============================================================
all_results_long <- bind_rows(lapply(names(all_results_by_scenario), function(scen) {
  df <- all_results_by_scenario[[scen]]
  df$scenario <- scen
  df
}))

# Save results
# save(all_results_long, file = here(res_dir, "gumbel_wholemax_res.Rdata"))
load(here(res_dir, "clust_gumbel_wholemax_res.Rdata"))
# ============================================================
# Summaries and quick plots
# ============================================================
# Convert to long format for plotting distributions of parameter estimates across R repetitions
all_long <- all_results_long %>%
  pivot_longer(cols = c(location, scale, shape), names_to = "parameter", values_to = "estimate")

# ============================================================
# Remove extreme outliers before plotting
# ============================================================
clean_long <- all_long

clean_long$keep <- TRUE   # logical flag

params <- unique(clean_long$parameter)

for (p in params) {
  idx <- clean_long$parameter == p
  est <- clean_long$estimate[idx]

  # compute quantile limits for *this* parameter
  q_low  <- quantile(est, 0.005, na.rm = TRUE)
  q_high <- quantile(est, 0.995, na.rm = TRUE)

  # mark as FALSE anything outside the limits
  clean_long$keep[idx] <- est >= q_low & est <= q_high
}

# keep filtered rows only
clean_long <- clean_long[clean_long$keep, ]

# remove the temporary flag
clean_long$keep <- NULL


# Example boxplot across scenarios (one point per repetition)
my_cols <- c("iid" = "#66c2a5", "theta_1.5" = "#fc8d62", "theta_2.5" = "#8da0cb")

ggplot(clean_long, aes(x = scenario, y = estimate, fill = scenario)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~parameter, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = my_cols) +
  labs(title = "Distribution of GEV parameter estimates over R repetitions (whole-sequence maxima)",
       subtitle = paste("Each repetition: generate", B, "maxima (each max from n =", n, "observations). R =", R),
       x = NULL, y = "Estimated parameter") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# You can also inspect numeric summaries:
all_results_long %>%
  group_by(scenario) %>%
  summarize(across(c(location, scale, shape), list(mean = mean, sd = sd), na.rm = TRUE))

load(here(plots_dir, "sim_plot.Rdata"))
print(p)
