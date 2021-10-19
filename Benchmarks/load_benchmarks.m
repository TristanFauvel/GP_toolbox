function [g, theta, model] = load_benchmarks(objective, details, benchmarks_table, rescaling)

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
    kernelname = char(benchmarks_table(benchmarks_table.fName == objective, :).kernelname);
    kernelfun = str2func(char(benchmarks_table(benchmarks_table.fName == objective, :).Kernel));
elseif ~isempty(details) && ~isempty(details.kernelname)
    kernelfun = str2func([details.kernelname, '_kernelfun']);
     kernelname = details.kernelname;
end

if isempty(details) || isempty(details.theta)
    theta = benchmarks_table(benchmarks_table.fName == objective, :).(kernelname);
    theta = theta{:};
elseif ~isempty(details) && ~isempty(details.theta)
    theta = details.theta;
end

     
model.hyp_lb = -10*ones(size(theta));
model.hyp_ub = 10*ones(size(theta));

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
 