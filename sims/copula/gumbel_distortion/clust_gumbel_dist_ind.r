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
n <- 100            # length of each sequence X_1...X_n
B <- 10             # how many maxima per repetition (M_1 ... M_B)
R <- 2              # how many repetitions of the whole experiment
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
save(all_results_long, file = here(res_dir, "clust_gumbel_wholemax_res.Rdata"))
