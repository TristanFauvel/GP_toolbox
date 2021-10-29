# GP_toolbox documentation
>Copyright Tristan Fauvel (2021)
>This software is distributed under the MIT License. Please refer to the file LICENCE.txt included for details

**GP_toolbox** is a **Matlab** toolbox for Gaussian process regression, classification and preference learning.The code is available at https://github.com/TristanFauvel/GP_toolbox. 
This toolbox was design so as to be highly flexible. The main reason why you may want to use this code instead of toolboxes such as the GPML toolbox is that this code is transparent and easy to modify. It also allows using state-of-the-art methods for approximate sampling from GP.

## Installation
  Simply add **GP_toolbox** to your Matlab path


## List of variables
- `xtrain` : design matrix
- `ytrain` : observations
- `xtest` : inputs for which a prediction is performed
- `post` : structure containing precomputed elements for inference (such as the inverse of the covariance matrix on the training set)
- theta : hyperparameters (theta.cov : covariance function hyperparameters, theta.mean : mean function hyperparameters).
- `regularization` : method used to regularize covariance matrices if needed, (either 'none' or 'nugget')
- `D`: dimension of the input space
- `meanfun` : prior mean function
- `kernelfun` : prior covariance function
- `ub`, `lb` : upper and lower bound of the space (assumed to be rectangular)
- `hyps` : structure containing informations about the hyperparameters (`ncov` : number of kernel hyperparameters, `nmean` : number of mean hyperparameters, `hyp_lb` and `hyp_ub` : bounds on the hyperparameters (in the form [bounds_on_kernel_hyps, bounds_on_mean_hyps]).

## GP models
A GP model is defined by the type of problem (regression, (binary) classification, preference), a kernel and a mean.
The fixed components of the model are stored within the object `model` corresponding to the class `gpmodel`.
There are three types of GP models implemented :
- GP regression: defined using the class `gp_regression_model`
- GP classification (for binary classification): defined using the class `gp_classification_model`
Both classes inherit from the class gpmodel.
- GP preference learning (based on pairwise comparisons) defined using the class `gp_preference_model`. This class inherits from `gp_classification_model`.


Any object from the class `gpmodel` has two main methods :
- `model_selection` : used to select the hyperparameters of the mean function and kernel.
- `prediction` : used for prediction.

For classification models, you need to specify:
- a link function (`normcdf()` or `logistic`).
- an approximation method (`'modeltype'`), either `'laplace`' (laplace approximation) or `'exp_prop'` (expectation propagation).

## Model selection

To perform model selection using maximum-likelihood estimation, you can use `model_selection(model, xtrain, ytrain, theta, update)`. update defines whether `‘all’` parameters should be learned, the kernel hyperparameters (`‘cov’`), or the mean (`‘mean’`) hyperparameters. The function returns the learned hyperparameters.


## Inference and  prediction
The method `prediction()` performs inference and prediction for inputs `xtest`.
Note that in the current implementation, the inputs are not automatically normalized, so it is better to pass min-max normalized inputs to this function.

### Prediction for regression :
- `mu_y` : latent mean
- `sigma2_y` : latent variance
- `Sigma2_y` : latent covariance
- `dmuy_dx`, `dsigma2_dx`, `dSigma2_dx`: derivatives with respect to the test input xtest.


### Prediction for classification and preference learning:

Observations `ctrain` are binary (0 or 1).
- `mu_c` : predictive class probability
-	`mu_y` : predictive mean
- `sigma2_y` :  predictive variance
- `Sigma2_y` : predictive covariance
- `var_muc` : variance of `link(f(x))` (a.k.a. epistemic uncertainty on c)  
- `dmuc_dx`, `dmuy_dx`, `dsigma2_dx`, `dSigma2_dx`, `dvar_muc_dx` : derivatives with respect to the test input xtest.




## Approximate sampling from GP

To perform approximate sampling, you can use the method `model.draw_sample_GP()`, by specifying the approximation type in the `approximation` structure.

The following methods are implemented to compute the features for finite dimensional kernel approximation:
 - [x] Sparse-spectrum approximation with random Fourier features (Lazaro-Gredilla et al, 2010) (`approximation.method = 'SSGP'`).
 - [x] Hilbert space method for reduced rank approximation (Solin & Särkkä, 2020) (`approximation.method = 'RRGP'`).

The following methods are implemented for approximate sampling:
  - [x] Weight-space approximation (Lazaro-Gredilla et al, 2010) (`approximation.decoupled_bases = 0`).
  - [x] Decoupled-bases approximation (Wilson et al, 2020) (`approximation.decoupled_bases = 1`).

In order to precompute the kernel approximation, use `model = approximate_kernel()` (note however, that if you change the model hyperparameters `theta`, you need to recompute the kernel approximation).


## Kernels

Various kernels are implemented and can be found in the /kernels subfolder. Functions computing the spectral density of stationary kernels are also implemented.

Additive kernels can be created using `sum_kernel()` and multiplied using `product_kernel()`.


## Note on preference learning

Preference learning with pairwise comparisons is essentially equivalent to binary classification. However, to account for the specific nature of the objective (preference) function, a specific preference kernel is used.
For a given base kernel `base_kernelfun`, the preference kernel is defined using `preference_kernelfun`:
`kernelfun= @(theta, xi, xj, training, reg)` `preference_kernelfun(theta, base_kernelfun, xi, xj, training, reg)`. In practice, you simply need to provide a base kernel and the corresponding preference kernel will be automatically computed.

 For improved uncertainty quantification, you can condition the value function on (x0, y0). To do so, you can define a kernel taking this conditioning into account using `conditioned_kernelfun`. In practice, you simply need to provide a `condition` when creating a preference learning model.

## General comments:
- The nugget regularization works by adding the smallest possible scalar to the covariance matrix diagonal. This is based on  Mohammadi et al, 2016. Note, however, that this can be slow for large matrices.
- When inference has already been performed for `(xtrain, ytrain)`, pass an inputs `post` to `prediction()` (otherwise, use`[]`). Precomputing `post`  can significantly accelerate `prediction()`.
- GP classification with non-zero mean is not  implemented.
- The **GP_toolbox** can be used with the **BO_toolbox** for Bayesian optimization.

## Reference
 If you use this software, please reference it as follows : Tristan Fauvel (2021) GP_toolbox, a Matlab toolbox for Gaussian process regression, classification and preference learning.

## Future features
 - [ ] Allow for non-zero mean in expectation-propagation approximation.
 - [ ] Allow for non-zero mean when using approximate sampling methods.
