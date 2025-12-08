functions {
  // --- HOFERT ET AL (2013) POLYNOMIAL LOGIC ---
  // (Same helper functions as before, kept for completeness)
  
  real log_abs_falling_factorial(real x, int n) {
    real res = 0;
    for (i in 0:(n-1)) res += log(abs(x - i));
    return res;
  }

  real sign_falling_factorial(real x, int n) {
    int neg_terms = 0;
    for (i in 0:(n-1)) {
      if ((x - i) < 0) neg_terms += 1;
    }
    return (neg_terms % 2 == 0) ? 1.0 : -1.0;
  }

  real log_sum_exp_signed(vector vals, vector signs) {
    int N = rows(vals);
    real max_val = max(vals);
    real sum_pos = 0;
    real sum_neg = 0;
    for (i in 1:N) {
      if (signs[i] > 0) sum_pos += exp(vals[i] - max_val);
      else sum_neg += exp(vals[i] - max_val);
    }
    real diff = sum_pos - sum_neg;
    if (diff <= 0) return negative_infinity(); // Catch numerical errors
    return max_val + log(diff);
  }

  real gumbel_high_dim_lpdf(vector u, real theta) {
    int d = rows(u);
    real alpha = 1.0 / theta;
    
    // Safety: Ensure u is strictly in (0,1)
    // We handle this in the model block, but this is a backup
    vector[d] safe_u;
    for(i in 1:d) safe_u[i] = fmin(fmax(u[i], 1e-15), 1.0 - 1e-15);

    vector[d] neg_log_u = -log(safe_u);
    vector[d] log_neg_log_u = log(neg_log_u);
    
    real t = sum(exp(theta * log_neg_log_u));
    real x = pow(t, alpha);
    
    vector[d] log_terms;
    vector[d] signs;
    
    for (j in 1:d) {
      real log_fall = log_abs_falling_factorial(alpha * j, d);
      real s_j_poly = sign_falling_factorial(alpha * j, d);
      real s_j_explicit = ( (d - j) % 2 == 0 ? 1.0 : -1.0 ) * s_j_poly;
      signs[j] = s_j_explicit;
      log_terms[j] = log_fall + j * log(x) + x - lgamma(j + 1) + poisson_lcdf(d - j | x);
    }
    
    real log_polynomial = log_sum_exp_signed(log_terms, signs);
    
    return d * log(theta) - x + (theta - 1) * sum(log_neg_log_u) 
           - d * log(t) - sum(log(safe_u)) + log_polynomial;
  }
}

data {
  int<lower=2> D;
  vector[D] Y;
}

parameters {
  real<lower=1.0> theta;
  real<lower=0, upper=1> w;
  ordered[2] mu;
  
  // FIX: Increased lower bound slightly to 0.1 to avoid the 
  // "cliff" where derivatives explode.
  vector<lower=0.1>[2] sigma; 
}

model {
  // --- Local Variable U (Not Saved) ---
  vector[D] U;

  // --- Priors ---
  theta ~ gamma(2, 1);
  w ~ beta(2, 2); // Slight preference for central weights
  mu ~ normal(0, 10);
  
  // Stronger regularization on sigma to prevent collapse
  sigma ~ gamma(4, 2); 

  // --- 1. Marginal Likelihood ---
  for (i in 1:D) {
    target += log_mix(w, 
                      normal_lpdf(Y[i] | mu[1], sigma[1]), 
                      normal_lpdf(Y[i] | mu[2], sigma[2]));
                      
    // --- 2. Compute U for Copula ---
    // We compute the CDF value. 
    real cdf_val = w * normal_cdf(Y[i]| mu[1], sigma[1]) 
                 + (1 - w) * normal_cdf(Y[i]| mu[2], sigma[2]);
    
    // SOFT CLAMPING (Nugget):
    // Instead of hard fmin/fmax, we mix with a tiny uniform 
    // to keep the gradient non-zero at the edges.
    // U_smooth = (1 - epsilon)*CDF + epsilon*0.5
    real eps = 1e-9;
    U[i] = (1.0 - 2.0*eps) * cdf_val + eps;
  }

  // --- 3. Copula Likelihood ---
  target += gumbel_high_dim_lpdf(U | theta);
}