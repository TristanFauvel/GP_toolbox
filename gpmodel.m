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
        xbest_norm = NaN;
        phi = []
        dphi_dx = [];
        ns = 0;
        context = false;
        grid
        grid_norm
        s0= 0.5;
        sdims
        xdims
        xgrid
        sgrid
        xgrid_norm
        sgrid_norm
        ongrid
    end
    methods
        function model = gpmodel(D, meanfun, kernelfun, regularization, hyps, lb, ub, kernelname, ns)
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
            model.ns = ns;
            model.sdims = 1:model.ns;
            model.xdims = (model.ns+1):model.D;

            if ns>0
                model.context = true;
            end
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
            if model.ongrid
                [~,  mu_y] =  model.prediction(theta, xtrain_norm, ctrain, model.grid_norm, post);
                mu_y = mean(mu_y, model.sdims);
                [ymax, id] =  max(mu_y);
                xmax = model.xgrid_norm(id);
            else
                opt_lb = model.lb_norm;
                opt_ub = model.ub_norm;
                opt_lb(model.sdims) = model.s0;
                opt_ub(model.sdims) = model.s0;
                init_guess = [];
                options.method = 'lbfgs';
                options.verbose = 1;
                ncandidates = 5;
                [xmax, ymax]  = optimize_AF(@(x)to_maximize_mean(model, theta, xtrain_norm, ctrain, x,post), opt_lb, opt_ub,ncandidates,init_guess, options, 'dims', model.xdims, 'objective', 'max');
                xmax = xmax(model.xdims);
            end

        end

        function model = approximate_kernel(model, theta, approximation)
            if strcmp(model.type, 'preference')
                [model.phi_pref, model.dphi_pref_dx, model.phi, model.dphi_dx] = sample_features_preference_GP(theta, model, approximation);
            else
                [model.phi, model.dphi_dx] = sample_features_GP(theta, model, approximation);

            end
        end

        function model = model_grid(model)
            grid_resolution = floor(5000^(1/model.D));
            model.grid_norm = myndgrid(linspace(0,1,grid_resolution), model.D)';
            model.grid = model.grid_norm.*(model.ub-model.lb)+model.lb;
        end

        function model = contextual_grid(model)
            grid_resolution = floor(5000^(1/(model.D-model.ns)));
            model.xgrid_norm = myndgrid(linspace(0,1,grid_resolution), model.D-model.ns)';

            grid_resolution = floor(5000^(1/(model.ns)));
            model.sgrid_norm = myndgrid(linspace(0,1,grid_resolution), model.ns)';

            model.xgrid = model.xgrid_norm.*(model.ub(model.xdims)-model.lb(model.xdims))+model.lb(model.xdims);
            model.sgrid = model.sgrid_norm.*(model.ub(model.sdims)-model.lb(model.sdims))+model.lb(model.sdims);
        end
    end
end
