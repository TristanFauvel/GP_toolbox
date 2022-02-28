function [g, theta, model] = load_benchmarks(objective, details, benchmarks_table, rescaling, type, varargin)
opts = namevaluepairtostruct(struct( ...
    'modeltype', 'exp_prop', ...
    'link', @normcdf...
    ), varargin);

UNPACK_STRUCT(opts, false)

obj = str2func(objective);

if isempty(details)
    obj = obj(rescaling);
else
    obj = obj(details.D, details.theta, details.kernelname, details.seed);
end
g = @(xx) obj.do_eval(xx);
xbounds = obj.xbounds;
D = obj.D;

if isempty(details) || isempty(details.kernelname)
    kernelname = char(benchmarks_table(benchmarks_table.fName == objective, :).Kernel_name);
    kernelfun = str2func(char(benchmarks_table(benchmarks_table.fName == objective, :).Kernel));
elseif ~isempty(details) && ~isempty(details.kernelname)
    kernelfun = str2func([details.kernelname, '_kernelfun']);
    kernelname = details.kernelname;
end

if isempty(details) || isempty(details.theta)
    h= benchmarks_table(benchmarks_table.fName == objective, :).(kernelname);
    theta.cov = h{:};
elseif ~isempty(details) && ~isempty(details.theta)
    theta.cov = details.theta;
end

meanfun = @constant_mean;
hyps.ncov_hyp =numel(theta.cov); % number of hyperparameters for the covariance function
hyps.nmean_hyp =1; % number of hyperparameters for the mean function
hyps.hyp_lb = -10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
hyps.hyp_ub = 10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
lb_norm = zeros(D,1);
ub_norm = ones(D,1);
lb = xbounds(:,1);
ub = xbounds(:,2);

% theta.mean = 0;
regularization = 'nugget';
ns = 0; %number of context variables
if strcmp(type, 'regression')
    model = gp_regression_model(D, meanfun, kernelfun, regularization, hyps, lb,ub, kernelname);
elseif strcmp(type, 'classification')
    model = gp_classification_model(D, meanfun, kernelfun, regularization, hyps, lb,ub, type, link, modeltype, kernelname, ns);
elseif strcmp(type, 'preference')
    condition.x0 = zeros(D,1);
    condition.y0 = 0;
    link = @normcdf;
    model = gp_preference_model(D, meanfun, kernelfun, regularization, hyps, lb,ub, type, link, modeltype, kernelname, condition, ns);
end

return
