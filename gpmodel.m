classdef gpmodel
    properties
        D
        meanfun
        kernelfun
        regularization
        hyp_lb
        hyp_ub
        ncov_hyp
        nmean_hyp
        ub
        lb
    end
    methods
        function model = gpmodel(D, meanfun, kernelfun, regularization, hyps, ub, lb)
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
        end
        
        function learned_hyp = model_selection(model, xtrain, ytrain, theta, update)
            switch update
                case 'cov'
                    model.hyp_lb((model.ncov_hyp+1):end)= theta.mean;
                    model.hyp_ub((model.ncov_hyp+1):end)= theta.mean;
                case 'mean'
                    model.hyp_lb(1:model.ncov_hyp) = theta.cov;
                    model.hyp_ub(1:model.ncov_hyp) = theta.cov;
            end
            %update_types = {'cov' , 'mean', 'all', 'none'};
            
            %else% ~any(strcmp(update_types, update))
            %    error('Specify a valid type of hyperparameters update')
            % end
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
    end
end
