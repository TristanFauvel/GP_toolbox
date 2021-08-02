function [g, theta, lb, ub, lb_norm, ub_norm, theta_lb, theta_ub, kernelfun, kernelname] = load_benchmarks(objective, kernelname, benchmarks_table, rescaling)
obj = str2func(objective);
obj = obj(rescaling);
g = @(xx) obj.do_eval(xx);
xbounds = obj.xbounds;
D = obj.D;

if isempty(kernelname)
    kernelname = char(benchmarks_table(benchmarks_table.fName == objective, :).Kernel_name);
    kernelfun = str2func(char(benchmarks_table(benchmarks_table.fName == objective, :).Kernel));
end

theta = benchmarks_table(benchmarks_table.fName == objective, :).(kernelname);
theta = theta{:};
theta_init = theta;
theta_lb = -10*ones(size(theta_init));
theta_ub = 10*ones(size(theta_init));

lb_norm = zeros(D,1);
ub_norm = ones(D,1);
lb = xbounds(:,1);
ub = xbounds(:,2);

return
% if strcmp(objective, 'forretal08')
%     xbounds = [0,1];
%     d=1; %dimension of the input space
%     if strcmp(kernelname, 'ARD')
%         theta = [2.8417; 7.7138]; %[3.0520, 0.6696]
%     end
%     xmax = 1;
%     scaling = 0.1;
%     bias =forretal08(xmax);
%     g = str2func(objective);
% elseif strcmp(objective, 'GP1d')
%     xbounds = [0,1];
%     d=1; %dimension of the input space
%     if strcmp(kernelname, 'ARD')
%         %         theta  = [4.0826,2]';
%         theta  = [3,2]'; %[5.8,2]';
%     end
%     scaling = 1;
%     rng(1)
%     xrange = linspace(xbounds(1),xbounds(2),100);
%     gen_kernelfun = @ARD_kernelfun;
%
%     gx =  mvnrnd(constant_mean(xrange,0), gen_kernelfun(theta, xrange,xrange)); %generate a function
%     gx = gx -min(gx);
%     g = sample_GP(theta, xrange, gx, kernelname,'RRGP', 1, 256, gen_kernelfun);
%
%
%     figure(); plot(gx)
% elseif strcmp(objective, 'grlee12')
%     xbounds = [0.5,2.5];
%     d=1; %dimension of the input space
%     if strcmp(kernelname, 'ARD')
%         theta  = [4.0826, 6.3708]';
%     end
%     scaling = 1;
%     %                         bias = grlee1(1);
%     g = str2func(objective);
%
% elseif strcmp(objective, 'levy')
%     d=2; %dimension of the input space
%     xbounds = [-10, 10; -10, 10];
%     if strcmp(kernelname, 'ARD')
%         theta = [3.9437,4.8886,13.3404]; %[4.9244,5.4415,6.8752]
%     end
%     scaling = 1;
%     %         bias = levy(1);
%     g = str2func(objective);
%
% elseif strcmp(objective, 'goldpr')
%     d=2; %dimension of the input space
%     xbounds = [-2, 2; -2, 2];
%     if strcmp(kernelname, 'ARD')
%         theta = [4.5784,7.5781, 15.0000]; %[4.3574, 6.5804, 10.0000]
%     end
%     scaling = 1;
%     %         bias = goldpr(1);
%     g = str2func(objective);
%
% elseif strcmp(objective, 'camel6')
%     d=2; %dimension of the input space
%     xbounds = [-3, 3; -2, 2];
%     if strcmp(kernelname, 'ARD')
%         theta = [4.0870, 2.3590, 8.6282];
%     end
%     scaling = 1;
%     %         bias = camel6(1);
%     g = str2func(objective);
% elseif strcmp(objective, 'Ursem_waves')
%     d= 2;
%      xbounds = [-0.9, 1.2; -1.2, 1.2];
%      if strcmp(kernelname, 'ARD')
%         theta = [4.7625, 4.1967,0.8135];
%      end
%     scaling = 1;
%     g = str2func(objective);
% end
