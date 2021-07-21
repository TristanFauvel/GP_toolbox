close all
add_gp_module

seed=1;
maxiter= 1000;

objectives = {'ackley'; 'beale'; 'boha1'; 'camel3'; 'camel6'; 'colville'; 'crossit'; ...
    'dixonpr'; 'drop'; 'egg'; 'forretal08'; 'goldpr'; 'griewank'; 'grlee12';  ...
    'hart3'; 'hart4'; 'hart6'; 'holder'; 'langer';'levy';'levy13';'perm0db';'permdb';...
    'powell'; 'rosen'; 'rothyp'; 'schaffer4'; 'schwef'; 'shekel'; 'shubert'; ...
    'spheref'; 'sumsqu'; 'trid'; 'Ursem_waves'};

kernelnames= {'ARD', 'Matern32', 'Matern52'};
kernelfuns= {'ARD_kernelfun', 'Matern32_kernelfun', 'Matern52_kernelfun'};

nk = numel(kernelnames);
N= numel(objectives);
load([pathname, '/Benchmarks/benchmarks_table.mat'], 'benchmarks_table')
T = benchmarks_table;
scores = zeros(N,nk);

theta.mean = 0;
        meanfun= @constant_mean;
for j = 1:N 
    bias = 0;
    objective = objectives{j};
    for i = 1:nk
        disp(['Function : ' , num2str(j), ', Kernel : ', num2str(i)]) 
        kernelname = kernelnames{i};
        kernelfun = str2func(kernelfuns{i});
        [g, theta.cov, lb, ub, lb_norm, ub_norm, theta_lb, theta_ub] = load_benchmarks(objective, kernelname, benchmarks_table);
        
        rng(seed)
        %% Initialize the experiment
        xtrain = rand_interval(lb,ub,'nsamples',maxiter);
        
        xtest = rand_interval(lb,ub,'nsamples',3*maxiter);

        ytrain = g(xtrain);
        ytrain_norm = ytrain - mean(ytrain);
        
        ytest = g(xtest);
        ytest = ytest - mean(ytest);

        %% Normalize data so that the bound of the search space are 0 and 1.
        xtrain_norm = (xtrain - lb)./(ub-lb);
        predictions = prediction(theta, xtrain_norm, ytrain_norm, xtest, kernelfun, meanfun);
        scores(j, i)  = mean((ytest(:) - predictions).^2);
    end
    [a,b] = min(scores(j,:),[], 2);
        best_kernel{j,1} = kernelfuns{b};
         best_kernel_name{j,1} = kernelnames{b};
end
benchmarks_table.Kernel = best_kernel;
benchmarks_table.Kernel_name = best_kernel_name;

save([pathname, '/Benchmarks/benchmarks_table.mat'], 'benchmarks_table')

% if strcmp(objective, 'forretal08')
%         theta = [4.7312 0.7403];
%
%     elseif strcmp(objective, 'grlee12')
%         theta  = [4.0826, 6.3708]';
%     elseif strcmp(objective, 'levy')
%         theta = [4.9244,5.4415,6.8752];
%     elseif strcmp(objective, 'goldpr')
%         theta = [4.3574, 6.5804, 10.0000];
%     elseif strcmp(objective, 'camel6')
%         theta = [4.0870, 2.3590, 8.6282];
%     end