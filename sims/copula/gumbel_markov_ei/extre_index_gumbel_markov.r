set.seed(42)

res_dir <- here("sims", "copula", "gumbel_markov_ei", "res")
plots_dir <- here(res_dir, "plots")

# -----------------------------
# Parameters
# -----------------------------
alpha <- 2.5
gumb <- gumbelCopula(alpha)

theta_theory <- 2^(1/alpha) - 1
cat("Theoretical extremal index:", theta_theory, "\n")

n <- 10000
burn <- 100
n_reps <- 50
thresholds <- seq(0.90, 0.997, length.out = 20)

theta_hat_reps <- matrix(NA, nrow = length(thresholds), ncol = n_reps)

# -----------------------------
# Tracking success/failure
# -----------------------------
attempts <- 0
successes <- 0

# -----------------------------
# Safe Markov chain simulation
# -----------------------------
rep_count <- 0

while (rep_count < n_reps) {

  attempts <- attempts + 1   # track attempts

  success <- TRUE
  U <- numeric(n + burn)
  U[1] <- runif(1)

  for (t in 1:(n + burn - 1)) {

    V <- runif(1)

    # --- safe conditional inversion ---
    next_val <- tryCatch(
      cCopula(cbind(U[t], V), copula = gumb, inverse = TRUE)[2],
      error = function(e) NA
    )

    if (is.na(next_val)) {
      success <- FALSE
      break
    }

    U[t + 1] <- next_val
  }

  # If chain failed, retry without counting as replicate
  if (!success) {
    cat("Replicate failed; retrying...\n")
    next
  }

  # Chain succeeded
  successes <- successes + 1

  rep_count <- rep_count + 1

  U <- U[(burn + 1):(n + burn)]   # drop burn-in

  theta_hat_reps[, rep_count] <- sapply(thresholds, function(u) {
    idx <- which(U[-length(U)] > u)
    if (length(idx) < 30) return(NA)
    1 - mean(U[idx + 1] > u)
  })

  cat("Replicate", rep_count, "completed successfully\n")
}

cat("All replicates complete.\n")

# -----------------------------
# Report success ratio
# -----------------------------
cat("\nTotal attempts:", attempts, "\n")
cat("Successful replicates:", successes, "\n")
cat("Success ratio:", round(successes / attempts, 3), "\n")

save(theta_hat_reps, file = here(res_dir, "theta_hat_reps.Rdata"))
# -----------------------------
# Average across replicates
# -----------------------------
theta_hat_mean <- rowMeans(theta_hat_reps, na.rm = TRUE)
theta_hat_sd   <- apply(theta_hat_reps, 1, sd, na.rm = TRUE)

# -----------------------------
# Plot results
# -----------------------------
plot(thresholds, theta_hat_mean, type = "b", pch = 16,
     ylim = c(0,1),
     xlab = "Threshold u", ylab = "Extremal Index θ(u)",
     main = sprintf("Gumbel Copula Markov Chain (alpha=%.2f) with %d replicates", alpha, n_reps))
arrows(thresholds, theta_hat_mean - theta_hat_sd,
       thresholds, theta_hat_mean + theta_hat_sd,
       angle = 90, code = 3, length = 0.02, col = "grey")

abline(h = theta_theory, col = "red", lwd = 2)
legend("bottomleft", legend = c("Empirical θ(u)", "Theoretical θ"),
       col = c("black","red"), pch = c(16, NA), lwd = c(1,2))
