%Copyright 2021 Tristan Fauvel
%This software is distributed under the MIT License. Please refer to the file LICENCE.txt included for details.
classdef gpmodel
    properties
        D
        meanfun
        kernelfun
        kernelname
        regularization
        hyp_lb
        hyp_ub
        ncov_hyp
        nmean_hyp
        ub
        lb
        ub_norm
        lb_norm
        max_muy = [];
     end
    methods
        function model = gpmodel(D, meanfun, kernelfun, regularization, hyps, lb, ub, kernelname)
            model.D = D;
            model.meanfun = meanfun;
            model.kernelfun = kernelfun;
            model.regularization = regularization;
            model.hyp_lb = hyps.hyp_lb;
            model.hyp_ub = hyps.hyp_ub;
            model.ncov_hyp = hyps.ncov_hyp;
            model.nmean_hyp = hyps.nmean_hyp;
            model.ub = ub;
            model.lb = lb;
            model.lb_norm = zeros(size(lb));
            model.ub_norm = ones(size(ub));
            model.kernelname = kernelname;
            
        end
        
        function learned_hyp = model_selection(model, xtrain, ytrain, theta, update)
            switch update
                case 'cov'
                    model.hyp_lb((model.ncov_hyp+1):end)= theta.mean;
                    model.hyp_ub((model.ncov_hyp+1):end)= theta.mean;
                case 'mean'
                    model.hyp_lb(1:model.ncov_hyp) = theta.cov;
                    model.hyp_ub(1:model.ncov_hyp) = theta.cov;
                case 'none'
                    learned_hyp = theta;
                    return
            end
            
            options_theta.method = 'lbfgs';
            hyp = multistart_minConf(@(hyp) model.minimize_negloglike(xtrain, ytrain, update, hyp), model.hyp_lb, model.hyp_ub, 10, [], options_theta);
            learned_hyp.cov = hyp(1:model.ncov_hyp);
            learned_hyp.mean = hyp(model.ncov_hyp+1:end);
        end
        
        %%
        function [negL, dnegL] =minimize_negloglike(model, xtrain, ytrain, update, hyp)
            hyp = hyp(:);
            if strcmp(update, 'none')
                negL = [];
                dnegL = zeros(1,model.ncov_hyp+model.nmean_hyp);
            else
                theta.cov = hyp(1:model.ncov_hyp);
                theta.mean = hyp((model.ncov_hyp+1):(model.ncov_hyp+model.nmean_hyp));
                [negL, dnegL] = model.negloglike(theta,  xtrain, ytrain);
                if strcmp(update, 'cov')
                    dnegL(model.ncov_hyp+1:end)=0;
                elseif strcmp(update, 'mean')
                    dnegL(1:model.ncov_hyp)=0;
                end
            end
        end
        
        function [xmax, ymax] =  maxmean(model, theta, xtrain_norm, ctrain, post)
            %% Return the maximum of the GP mean
            init_guess = model.max_muy;
            options.method = 'lbfgs';
            options.verbose = 1;
            ncandidates = 5;
            [xmax, ymax] = multistart_minConf(@(x) model.to_maximize_mean(theta, xtrain_norm, ctrain, x, post), ...
                model.lb_norm,  model.ub_norm, ncandidates, init_guess, options, 'objective', 'max');
        end
    end
end
