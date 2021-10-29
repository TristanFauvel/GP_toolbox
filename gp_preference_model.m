classdef gp_preference_model < gp_classification_model
    properties
        base_kernelfun
        condition
    end
    methods
        function model = gp_preference_model(D, meanfun, base_kernelfun, regularization, hyps, lb, ub, type, link, modeltype, kernelname, condition)
            if isempty(condition)
                kernelfun = @(theta, xi, xj, training, reg) preference_kernelfun(theta, base_kernelfun, xi, xj, training, reg);
            else
                kernelfun= @(theta, xi, xj, training, reg) conditional_preference_kernelfun(theta, base_kernelfun, xi, xj, training, reg, condition.x0);
            end

            model = model@gp_classification_model(D, meanfun, kernelfun, regularization, hyps, lb, ub, type, link, modeltype, kernelname);
            model.base_kernelfun = base_kernelfun;
            model.condition = condition;
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
    end
end

