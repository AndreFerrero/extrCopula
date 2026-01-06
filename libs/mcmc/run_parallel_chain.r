run_parallel_chains <- function(
  log_target, # log-posterior function
  init_values, # list of initial points (one per chain)
  n_iter,
  proposal,
  n_cores,
  adapt = NULL,
  transform = NULL, # inverse transform (e.g. g_inv_dep)
  export = NULL # objects/functions to export to workers
) {
    library(parallel)
    library(coda)

    # Create cluster
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)

    # Load required packages on workers
    clusterEvalQ(cl, {
        library(coda)
    })

    # Export needed objects to workers
    clusterExport(
        cl,
        varlist = unique(c(
            "run_chain",
            "n_iter",
            "log_target",
            "proposal",
            "adapt",
            "transform",
            "param_map",
            export
        )),
        envir = environment()
    )

    # Function run on each worker
    run_single_chain <- function(init) {
        res <- run_chain(
            log_target = log_target,
            init       = init,
            n_iter     = n_iter,
            proposal   = proposal,
            adapt      = adapt
        )

        samples <- if (!is.null(transform)) {
            t(apply(res$samples, 1, transform))
        } else {
            res$samples
        }

        mcmc(samples)
    }

    # Run chains in parallel
    chains <- parLapply(cl, init_values, run_single_chain)

    # Return mcmc.list
    mcmc.list(chains)
}
