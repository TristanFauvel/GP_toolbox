function [g_mu_y,  dmuy_dx] = to_maximize_mean_bin_GP(theta, xtrain_norm, ctrain, x, kernelfun, modeltype, post)

if any(isnan(x(:)))
    error('x is NaN')
end
[d,n]= size(x);

regularization = 'nugget';

if isempty(post)
    warning('Precomputing the approximate posterior is more efficient')
end
[~,  g_mu_y, ~, ~, ~, dmuy_dx] = prediction_bin(theta, xtrain_norm, ctrain, x, kernelfun, modeltype, post, regularization);


g_mu_y = -g_mu_y;
% dmuy_dx= -squeeze(dmuy_dx(:,:,1:d));
dmuy_dx= -squeeze(dmuy_dx);
% dmuy_dx = dmuy_dx(2:end);
return

n= 100;
x = linspace(0,1,n);
s0 = xtrain_norm(1,end);
[mu_c,  g_mu_y]= prediction_bin(theta, xtrain_norm, ctrain, [s0*ones(1,n);x], kernelfun, modeltype, post, regularization);

figure()
plot(x, g_mu_y)
figure()
plot(x, mu_c)