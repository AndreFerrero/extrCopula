# Code Library

This folder contains reusable code with a modular structure, with the intention of providing faster and easier implementation of new features.

`packages.r` is sourced whenever a function requires specific packages.

## MCMC

Here we keep all the MCMC related code.

- `adaptation`: Adaptation schemes
- `engines`: MCMC core step structure (e.g. Metropolis Hastings step)
- `proposals`: Proposal distributions (e.g. Multivariate gaussian)
- `run_chain.r`: Single chain loop combining previous pieces

## Models

Here we store model objects and useful builders

- `builders`: Builder functions that use model objects functionalities (Simulator with copula + margin choice, Log-posterior to be used in MCMC,...)
- `copulas`: Copula objects with helpers functions related to the copula choice (simulate_u, log_density, log_prior...)
- `margins`: Margin objects with related functionalities
