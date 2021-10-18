function [g_mu_y,  dmuy_dx] = to_maximize_mean_bin_GP(theta, xtrain_norm, ctrain, x,model, post)

if any(isnan(x(:)))
    error('x is NaN')
end
 
 
if isempty(post)
    warning('Precomputing the approximate posterior is more efficient')
end
[~,  g_mu_y, ~, ~, ~, dmuy_dx] = model.prediction(theta, xtrain_norm, ctrain, x, post);


g_mu_y = -g_mu_y;
dmuy_dx= -squeeze(dmuy_dx);

return

