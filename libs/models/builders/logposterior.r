source("libs/packages.R")

#' Joint posterior with full likelihood
#' 
#' @param copula Copula object
#' @param margin Margin object
#' @param param_map named list (margin and copula)containing parameter references for each of the 2 object
#' @param inverse_transform function to transform parameters in the original space
#' @param log_jacobian jacobian function of the transformation
#' 
#' @return Log posterior density (up to normalising constant)
#' 
build_logposterior <- function(copula, margin, param_map, data,
                               inverse_transform = NULL,
                               log_jacobian = NULL,
                               margin_prior = NULL,
                               copula_prior = NULL) {
  
  function(param_init) {
    
    # Apply inverse transform if given
    # param is always the parameter in the original space
    # no transformation leaves param = param_init
    param <- if (!is.null(inverse_transform)) inverse_transform(param_init) else param_init
    
    param_m <- param[param_map$margin]
    param_c <- param[param_map$copula]
    
    # --- Prior ---
    logprior_m <- do.call(margin$log_prior, c(list(param_m), margin_prior))
    logprior_c <- do.call(copula$log_prior, c(list(param_c), copula_prior))
    logprior <- logprior_m + logprior_c
    
    if (!is.finite(logprior)) return(-Inf)
    
    # 2. Likelihood (copula + marginal)
    u <- margin$cdf(data, param_m)
    u <- pmin(pmax(u, 1e-8), 1 - 1e-8)  # numerical stability
    
    loglik <- sum(margin$log_density(data, param_m)) +
          copula$log_density(u, param_c)
    
    # 3. Jacobian adjustment if in transformed space
    logjac <- if (!is.null(log_jacobian)) log_jacobian(param_init) else 0
    
    logprior + loglik + logjac
  }
}
