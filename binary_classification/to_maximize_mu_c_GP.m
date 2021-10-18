function [mu_c,  dmuc_dx] = to_maximize_muc_GP(theta, xtrain_norm, ctrain, x,model, post)

if any(isnan(x(:)))
    error('x is NaN')
end


if isempty(post)
    warning('Precomputing the approximate posterior is more efficient')
end
[mu_c,  mu_y, sigma2_y, Sigma2_y, dmuc_dx] = model.prediction(theta, xtrain_norm, ctrain, x, post);


mu_c= -mu_c;
dmuc_dx= -squeeze(dmuc_dx);

return

