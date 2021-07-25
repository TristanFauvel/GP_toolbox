function [g_mu_y,  dmuy_dx] = to_maximize_value_function(theta, xtrain_norm, ctrain, x, kernelfun,x0, modeltype, post)

if any(isnan(x(:)))
    error('x is NaN')
end
[D,n]= size(x);

regularization = 'nugget';
[~,  g_mu_y, ~, ~, ~, dmuy_dx] = prediction_bin(theta, xtrain_norm, ctrain, [x;x0*ones(1,n)], kernelfun, modeltype, post, regularization);


g_mu_y = -g_mu_y;
dmuy_dx= -squeeze(dmuy_dx(1:D,:,:));

return