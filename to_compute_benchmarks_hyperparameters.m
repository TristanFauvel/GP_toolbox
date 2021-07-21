clear all
close all
add_modules
pathname = 'C:\Users\tfauvel\Documents\Retinal_prosthetic_optimization\Retinal_Prosthetic_Optimization_Code\Refractive_error_measurement';
cd(pathname)
addpath(genpath(pathname))
figures_folder = pathname;
graphics_style_paper;

seed=1;
maxiter= 500 ;
rng(seed)

n=80;
sampling_method = 'SSGP';

options_theta.method = 'lbfgs';
options_theta.verbose = 1;

update_period = 15000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ninit = 15; % number of time steps before starting using the hyperparameters
nrepets = 1; %20;
seeds=1:nrepets;

gmax = NaN(1,nrepets);
xmax = NaN(1,nrepets);
show_plots = 0;

save_data = 1;
objectives = {'forretal08', 'grlee12', 'levy', 'goldpr', 'camel6'};

for j = 1:5
    bias = 0;
    objective = objectives{j};
    if strcmp(objective, 'forretal08')
        xbounds = [0,1];
        d=1; %dimension of the input space
        theta = [4.7312 0.7403];
        xmax = 1;
        scaling = 0.1;
        bias =forretal08(xmax);
        
    elseif strcmp(objective, 'grlee12')
        xbounds = [0.5,2.5];
        d=1; %dimension of the input space
        theta  = [4.0826, 6.3708]';
        scaling = 1;
        %                         bias = grlee1(1);
    elseif strcmp(objective, 'levy')
        d=2; %dimension of the input space
        xbounds = [-10, 10; -10, 10];
        theta = [4.9244,5.4415,6.8752];
        scaling = 1;
        %         bias = levy(1);
    elseif strcmp(objective, 'goldpr')
        d=2; %dimension of the input space
        xbounds = [-2, 2; -2, 2];
        theta = [4.3574, 6.5804, 10.0000];
        scaling = 1;
        %         bias = goldpr(1);
    elseif strcmp(objective, 'camel6')
        d=2; %dimension of the input space
        xbounds = [-3, 3; -2, 2];
        theta = [4.0870, 2.3590, 8.6282];
        scaling = 1;
        %         bias = camel6(1);
    end
    
    theta(end) = theta(end) +2*log(scaling);
    theta_init = theta;
    theta_lb = -10*ones(size(theta_init));
    theta_ub = 10*ones(size(theta_init));
    
    nd = d+1;
    lb_norm = zeros(nd,1);
    ub_norm = ones(nd,1);
    
    
    
    g = str2func(objective);
    g = @(x) scaling*(g(x) -bias);
    
    sbound = [0,1];
    
    
    lb = [sbound(1);xbounds(:,1)];
    ub = [sbound(2);xbounds(:,2)];
    
    kernelfun = @Refractive_kernelfun;
    
    nalt = 2;
    link = @normcdf; %inverse link function for the classification model
    b = norminv(1/nalt);
    f = @(s,x) normcdf(g(x).*s + b);
    modeltype = 'exp_prop'; % Approximation model
    
    kernelname = 'ARD';
    nacq = 2; % number of time steps before starting using the acquisition function
    
    theta= theta_init;
    
    
    rng(seed)
    %% Initialize the experiment
    xtrain = rand_interval(lb,ub,'nsamples',maxiter);
    %Generate a binary sample
    ctrain = f(xtrain(1,:), xtrain(2,:))>rand(1,maxiter);
    
    %% Normalize data so that the bound of the search space are 0 and 1.
    xtrain_norm = (xtrain - lb)./(ub-lb);
    %Compute the maximum of the value function according to the model
    
    options=[];
    init_guess = theta;
    theta = multistart_minConf(@(hyp)negloglike_bin(hyp, xtrain_norm, ctrain, kernelfun, 'modeltype', modeltype), theta_lb, theta_ub,20, init_guess, options_theta);
    disp([objective, ': '])
    disp(theta)
    
    [mu_c, mu_y, sigma2_y] = prediction_bin(theta, xtrain, ctrain, inputs_test, kernelfun, 'modeltype', modeltype);
    fig=figure('units','centimeters','outerposition',1+[0 0 16 16]);
    fig.Color =  [1 1 1];
    x = xtrain;
    c = ctrain;
    subplot(1,3,1)
    imagesc(s_range,xrange,reshape(mu_c,n,n)'); hold on;
    scatter(x(1,c == 0), x(2,c == 0),markersize, 'k', 'filled'); hold on;
    scatter(x(1,c == 1), x(2,c == 1),markersize,  'r', 'filled'); hold off;
    xlabel('Size')
    ylabel('Correction')
    set(gca,'YDir','normal')
    colorbar
    box off
    pbaspect([1 1 1])
    
    subplot(1,3,2)
    imagesc(s_range,xrange,reshape(mu_y,n,n)'); hold on;
    scatter(x(1,c == 0), x(2,c == 0),markersize, 'k', 'filled'); hold on;
    scatter(x(1,c == 1), x(2,c == 1),markersize,  'r', 'filled'); hold off;
    xlabel('Size')
    ylabel('Correction')
    set(gca,'YDir','normal')
    set(gca,'YDir','normal')
    colorbar
    box off
    pbaspect([1 1 1])
    subplot(1,3,3)
    imagesc(s_range,xrange,reshape(sigma2_y,n,n)'); hold on;
    scatter(x(1,c == 0), x(2,c == 0),markersize, 'k', 'filled'); hold on;
    scatter(x(1,c == 1), x(2,c == 1),markersize,  'r', 'filled'); hold off;
    xlabel('Size')
    ylabel('Correction')
    set(gca,'YDir','normal')
    set(gca,'YDir','normal')
    colorbar
    box off
    pbaspect([1 1 1])

end