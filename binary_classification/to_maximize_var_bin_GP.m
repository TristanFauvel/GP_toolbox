function [var_muc, dvar_muc_dx] = to_maximize_var_bin_GP(theta, xtrain_norm, ctrain, x,model, post)

if any(isnan(x(:)))
    error('x is NaN')
end
[d,n]= size(x);

regularization = 'nugget';

if isempty(post)
    warning('Precomputing the approximate posterior is more efficient')
end
[output1,  mu_y, sigma2_y, Sigma2_y, dmuc_dx, dmuy_dx, dsigma2y_dx, dSigma2y_dx, var_muc, dvar_muc_dx,post] = ...
    model.prediction(theta, xtrain_norm, ctrain, x, post);


var_muc = -var_muc;
dvar_muc_dx= -squeeze(dvar_muc_dx);
return
