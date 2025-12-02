handy_dir <- here("sims", "common")
source(here(handy_dir, "handy_funs.r"))
# set.seed(46)

n <- 10000
alpha <- 2
theta <- 1.5

Gcop <- gumbelCopula(param = theta, dim = n)

X_v <- qfrechet(rGumbV(n, theta), shape = alpha)
X_c <- rCopFrechet(alpha, Gcop)

# par(mfrow = c(1,2))
# plot(density(X_v))
# plot(density(X_c))
# par(mfrow = c(1,1))

# Compute densities of your simulated data
dens_v <- density(X_v)
dens_c <- density(X_c)

# Plot the first density
plot(dens_v, col = "blue", lwd = 2, main = "Gumbel - Frechet Samplers", xlab = "Value", ylab = "Density")

# Add the second density
lines(dens_c, col = "red", lwd = 2)

# Overlay the theoretical Frechet density
x_vals <- seq(min(c(X_v, X_c)), max(c(X_v, X_c)), length.out = 1000)
dens_theoretical <- dfrechet(x_vals, shape = alpha)  # Assuming dfrechet exists from your handy_funs
lines(x_vals, dens_theoretical, col = "darkgreen", lwd = 2, lty = 2)

# Add legend
legend("topright", 
       legend = c("V sampler", "Copula package", "Frechet(2)"),
       col = c("blue", "red", "darkgreen"), 
       lwd = 2, lty = c(1,1,2))


gum_psi <- function(t, theta){
    exp(- t ^ (1/theta))
}

# set.seed(46)
Gcop <- gumbelCopula(param = theta, dim = n)

V <- rstable(
    n = 1,
    alpha = 1/theta,
    beta = 1,
    gamma = cospi(1/(2 * theta))^theta,
    delta = 0,
    pm = 1
)

E <- rexp(n, 1) 

psi_E <- gum_psi(E/V, theta) 


par(mfrow = c(1,2))
hist(psi_E)
hist(rCopula(1, Gcop))
par(mfrow = c(1,1))

