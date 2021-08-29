# GP_toolbox
 A Matlab toolbox for Gaussian process regression, classification and preference learning

## Reference
 If you use this software, please reference it as follows : Tristan Fauvel (2021) GP_toolbox, a Matlab toolbox for Gaussian process regression, classification and preference learning.

## Features
 * Gaussian process regression
 * Gaussian process classification with Laplace approximation or Expectation-Propagation
 * Gaussian process preference learning using conditional preference kernel
 * Approximate sampling from posterior GP :
   * Finite dimensional stationary kernel approximation using Sparse-Spectrum approximation (Lazaro-Gredilla et al, 2010)
   * Finite dimensional stationary kernel approximation using Hilbert-space methods (Solin et al, 2010)
   * Weight-space approximate sampling (Lazaro-Gredilla et al, 2010) or decoupled-bases approximate sampling (Wilson et al, 2020)

## Installation
   * Simply add GP_toolbox to your Matlab path

## User guide
   * Each /Examples subfolder includes detailed description about how to use the toolbox. 
   * You can use GP_toolbox with BO_toolbox for Bayesian optimization (https://github.com/TristanFauvel/BO_toolbox)


## License
   This software is distributed under the MIT License. Please refer to the file LICENCE.txt included for details.

## Future features
- [ ] Allow for non-zero mean in expectation-propagation approximation.
- [ ] Allow for non-zero mean when using approximate sampling methods.
