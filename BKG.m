function [new_x, new_x_norm] = BKG(theta, xtrain_norm, ctrain,model, post, approximation)
if ~strcmp(model.modeltype, 'laplace')
    error('This acquisition function is only implemented with Laplace approximation')
end
init_guess = [];
options.method = 'lbfgs';
options.verbose = 1;
ncandidates = 10;
[xbest, ybest] = multistart_minConf(@(x)to_maximize_mean_bin_GP(theta, xtrain_norm, ctrain, x, model, post), model.lb_norm, model.ub_norm, ncandidates, init_guess, options);
ybest = -ybest;

c0 = [ctrain, 0];
c1 = [ctrain,1];
regularization = 'nugget';
new_x_norm  = multistart_minConf(@(x)knowledge_grad(theta, xtrain_norm, ctrain, x,model, post, c0, c1, xbest, ybest, model.lb_norm, model.ub_norm), model.lb_norm, model.ub_norm, ncandidates, init_guess, options);

new_x = new_x_norm.*(model.ub-model.lb) + model.lb;
end
 