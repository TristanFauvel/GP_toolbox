classdef gp_preference_model < gp_classification_model
    properties
        condition
         phi_pref= []
        dphi_pref_dx = [];
        base_kernelfun = [];
    end
    methods
        function model = gp_preference_model(D, meanfun, base_kernelfun, regularization, hyps, lb, ub, type, link, modeltype, kernelname, condition, ns)
            if isempty(condition)
                kernelfun = @(theta, xi, xj, training, reg) preference_kernelfun(theta, base_kernelfun, xi, xj, training, reg);
            else
                kernelfun= @(theta, xi, xj, training, reg) conditional_preference_kernelfun(theta, base_kernelfun, xi, xj, training, reg, condition.x0);
            end

            model = model@gp_classification_model(D, meanfun, kernelfun, regularization, hyps, lb, ub, type, link, modeltype, kernelname, ns);
            model.condition = condition;
            model.base_kernelfun = base_kernelfun;
        end

        function [sample_f, dsample_f_dx, decomposition] = draw_sample_GP(model, theta, xtrain, ctrain, approximation, post)
            % The decomposition corresponds to the prior and update terms in
            % the decoupled bases approximation.
            if isempty(model.phi)
                model = approximate_kernel(model, theta, approximation);
            end

            approximation.phi = model.phi;
            approximation.dphi_dx = model.dphi_dx;

            approximation.phi_pref = @(x)  approximation.phi(x(1:model.D,:)) -  approximation.phi(x(model.D+1:end,:));
            approximation.dphi_pref_dx = @(x)  approximation.dphi_dx(x(1:model.D,:));


            [sample_f, dsample_f_dx, decomposition] = sample_binary_GP_precomputed_features(xtrain, ctrain, theta,model, approximation, post);
        end

        function [sample_normalized, sample] = sample_max_GP(model, approximation, xtrain_norm, ctrain, theta, post)
            options.method = 'lbfgs';
            options.verbose = 1;
            ncandidates = 5;
            approximation.phi = model.phi;
            approximation.dphi_dx = model.dphi_dx;
            approximation.phi_pref = model.phi_pref;
            approximation.dphi_pref_dx = model.dphi_pref_dx;
            [sample_g, dsample_g_dx] = sample_value_GP_precomputed_features(approximation, theta, xtrain_norm, ctrain, model, post);

            init_guess = [];
            new = multistart_minConf(@(x)deriv(x,sample_g, dsample_g_dx), model.lb_norm, model.ub_norm, ncandidates,init_guess, options, 'objective', 'max');

            %Careful here: the goal is to maximize the value function (the problem
            %is a maximization problem): deriv takes the opposite of sample_g
            sample = new.*(model.ub -model.lb) + model.lb;
            sample_normalized= new;
        end
    end
end

