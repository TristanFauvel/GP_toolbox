function [prior_mean, dprior_mean_dtheta, dprior_mean_dx] =constant_mean(x,hyp)
if numel(hyp)>1
error('The number of hyperparameters for this mean function is 1');
end
prior_mean = hyp*ones([1,size(x,2)])';
dprior_mean_dx =zeros(size(x))';
dprior_mean_dtheta = ones([1,size(x,2)]);
return