function [g_mu_y,  dmuy_dx] = to_maximize_value_function(theta, xtrain_norm, ctrain, x, model, post)

if any(isnan(x(:)))
    error('x is NaN')
end
[D,n]= size(x);

[~,  g_mu_y, ~, ~, ~, dmuy_dx] = prediction_bin(theta, xtrain_norm, ctrain, [x; model.condition.x0*ones(1,n)], model, post);


g_mu_y = -g_mu_y;
dmuy_dx= -squeeze(dmuy_dx(1:D,:,:));

return