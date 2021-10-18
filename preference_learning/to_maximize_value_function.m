%Copyright 2021 Tristan Fauvel
%This software is distributed under the MIT License. Please refer to the file LICENCE.txt included for details.

function [g_mu_y,  dmuy_dx] = to_maximize_value_function(theta, xtrain_norm, ctrain, x, model, post)
% Used in conjunction with minConf to maximize the latent value function of the preference GP model. 
if any(isnan(x(:)))
    error('x is NaN')
end
[D,n]= size(x);

[~,  g_mu_y, ~, ~, ~, dmuy_dx] = model.prediction(theta, xtrain_norm, ctrain, [x; model.condition.x0*ones(1,n)], post);


g_mu_y = -g_mu_y;
dmuy_dx= -squeeze(dmuy_dx(1:D,:,:));

return
