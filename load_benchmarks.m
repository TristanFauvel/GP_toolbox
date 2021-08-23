function [g, theta, model] = load_benchmarks(objective, kernelname, benchmarks_table, rescaling)
obj = str2func(objective);
obj = obj(rescaling);
g = @(xx) obj.do_eval(xx);
xbounds = obj.xbounds;
D = obj.D;

if isempty(kernelname)
    kernelname = char(benchmarks_table(benchmarks_table.fName == objective, :).Kernel_name);
    kernelfun = str2func(char(benchmarks_table(benchmarks_table.fName == objective, :).Kernel));
else
    kernelfun = str2func([kernelname, '_kernelfun']);
end

theta = benchmarks_table(benchmarks_table.fName == objective, :).(kernelname);
theta = theta{:};
theta_init = theta;
model.theta_lb = -10*ones(size(theta_init));
model.theta_ub = 10*ones(size(theta_init));

model.lb_norm = zeros(D,1);
model.ub_norm = ones(D,1);
model.lb = xbounds(:,1);
model.ub = xbounds(:,2);

model.kernelname = kernelname; 
model.kernelfun = kernelfun;
model.meanfun = @constant_mean;
model.regularization = 'nugget';
model.D = D;
return
 