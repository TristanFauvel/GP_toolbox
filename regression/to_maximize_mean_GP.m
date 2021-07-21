function [mu_y,  dmuy_dx] = to_maximize_mean_GP(theta, xtrain, ytrain, x, kernelfun, post)

if any(isnan(x(:)))
    error('x is NaN')
end

[mu_y, sigma2_y, dmu_dx] =  prediction(theta, xtrain, ytrain, x, kernelfun, meanfun, 'post', post);


mu_y = -mu_y;
dmuy_dx= -squeeze(dmu_dx);

return

