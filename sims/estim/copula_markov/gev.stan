functions {
  real gev_lpdf(real x, real mu, real sigma, real xi) {
    real z;
    if (xi != 0) {
      z = 1 + xi * (x - mu) / sigma;
      if (z <= 0) return negative_infinity(); // support check
      return -log(sigma) - (1/xi + 1) * log(z) - pow(z, -1/xi);
    } else {
      // Gumbel limit xi -> 0
      return -log(sigma) - (x - mu)/sigma - exp(-(x - mu)/sigma);
    }
  }
}
data {
  int<lower=1> N;       // number of block maxima
  array[N] real M;            // block maxima
}
parameters {
  real mu;
  real<lower=0> sigma;
  real xi;
}
model {
  mu ~ normal(0, 10);
  sigma ~ normal(0, 10);
  xi ~ normal(0, 5);

  for (n in 1:N)
    target += gev_lpdf(M[n] | mu, sigma, xi);
}
