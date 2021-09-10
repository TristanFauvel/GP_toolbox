function [g,  dgdx] = to_maximize_averaged_mean_bin_GP(theta, xtrain_norm, ctrain, x, model, post)
% Average over the contexts s
if any(isnan(x(:)))
    error('x is NaN')
end
 
 
if isempty(post)
    warning('Precomputing the approximate posterior is more efficient')
end

lb = model.lb_norm(1:model.ns);
ub = model.ub_norm(1:model.ns);

fun = @(s) mean_binGP(theta, xtrain_norm, ctrain, [s;x*ones(1,size(s,2))], model, post);
g = -integral(fun,lb,ub, 'ArrayValued', true,'RelTol',1e-4);

fun = @(s) dmean_binGP(theta, xtrain_norm, ctrain, [s;x*ones(1,size(s,2))], model, post);
dgdx = -integral(fun,lb,ub, 'ArrayValued', true,'RelTol', 1e-4);

 end

function g_mu_y = mean_binGP(theta, xtrain_norm, ctrain, x, model, post)

[~,  g_mu_y] = prediction_bin(theta, xtrain_norm, ctrain, x, model, post);

end
function  dmuy_dx= dmean_binGP(theta, xtrain_norm, ctrain, x, model, post)

[~,  ~, ~, ~, ~, dmuy_dx] = prediction_bin(theta, xtrain_norm, ctrain, x, model, post);

dmuy_dx = squeeze(dmuy_dx);
dmuy_dx = dmuy_dx((model.ns+1):end);
end