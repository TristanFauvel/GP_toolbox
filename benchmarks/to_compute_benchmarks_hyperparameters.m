close all
add_gp_module

seed=1;
maxiter= 1000;

options_theta.method = 'lbfgs';
options_theta.verbose = 1;

rescaling = 1;
binary = 1;

if rescaling == 1
    if binary == 1
        filename = [pathname, '/Benchmarks/benchmarks_table_rescaled_bin.mat'];
    else
        filename = [pathname, '/Benchmarks/benchmarks_table_rescaled.mat'];
    end
else
    if binary == 1
        filename = [pathname, '/Benchmarks/benchmarks_table_bin.mat'];
    else
        filename = [pathname, '/Benchmarks/benchmarks_table.mat'];
    end
end


load(filename, 'benchmarks_table')
objectives = benchmarks_table.fName;
kernelnames= {'ARD', 'Matern32', 'Matern52'};
kernelfuns= {'ARD_kernelfun', 'Matern32_kernelfun', 'Matern52_kernelfun'};

nk = numel(kernelnames);
N= numel(objectives);
T = benchmarks_table;
for j = 1:N %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bias = 0;
    objective = char(objectives(j));
    for i = 1:nk
        disp(['Function : ' , num2str(j), ', Kernel : ', num2str(i)])
        kernelname = kernelnames{i};
        kernelfun = str2func(kernelfuns{i});
        [g, theta_init, model] = load_benchmarks(objective, kernelname, benchmarks_table, rescaling);
        theta= theta_init;
        
        rng(seed)
        %% Initialize the experiment
        xtrain = rand_interval(model.lb ,model.ub ,'nsamples',maxiter);
        ytrain = g(xtrain);
        if size(ytrain, 2) ~= maxiter
            error('ouput_size')
        end
        %% Normalize data so that the bound of the search space are 0 and 1.
        xtrain_norm = (xtrain -model.lb)./(model.ub - model.lb);
        ytrain_norm = ytrain - mean(ytrain);
        
        %Compute the maximum of the value function according to the model
        
        init_guess = [];
        ncov_hyp = numel(model.hyp_lb);
        
        hyp.cov = 10*ones(ncov_hyp,1);%rand(ncov_hyp,1);
        hyp.mean = 0 ;
        
        
        update = 'cov';
        meanfun= @constant_mean;
        nmean_hyp = 1;
        theta = multistart_minConf(@(hyp)minimize_negloglike(hyp, xtrain_norm, ytrain_norm, kernelfun, meanfun, ncov_hyp, nmean_hyp, update), [model.hyp_lb; 0], [model.hyp_ub; 0],20, init_guess, options_theta);
        theta = theta(1:end-1);
        
        %         if strcmp(kernelname,'ARD')
        T(T.fName == objectives(j), i + 3) = mat2cell(theta, size(theta,1),size(theta,2));
        %         else
        %             T(T.fName == objectives{i}, i + 3) = theta(:)';
        %         end
    end
end
benchmarks_table = T;
% save(filename, 'benchmarks_table')
