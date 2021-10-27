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
        
    end
end

 