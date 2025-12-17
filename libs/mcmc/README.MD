# MCMC Scheme

The logic of this folder is to have one script running the chain loop (`run_chain.r`), with the inner functionalities defined outside for modularity.

## Chain logic

- `log_target`: Defines the log posterior of the model. This should be built using builders in `models` folder. Approximations of the posterior using Bayesian Synthetic Likelihood is allowed by defining the appropriate approximation scheme

- `init`: Initial parameters of the chains. If a transformation is desired, these are the parameters in the transformed space. Jacobian transformation should be implemented in the `log_target` function.

- `proposal`: Proposal distribution at each step. Usually this is a MVN.

- `engine_step`: Step structure, default is Metropolis Hastings.

- `adapt`: Adaptation scheme to be used.

To allow for adaptation, the chain uses an internal variable `state` which contains information on the proposal distribution. Different adaptation schemes require and update different information, and should enforce the creation of necessary components at the first iteration.