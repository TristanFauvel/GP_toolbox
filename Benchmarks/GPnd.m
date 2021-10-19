classdef GPnd
    properties
        xbounds
        D
        name = 'GP';
        b = 0.5;
        opt = 'max';
        theta
        kernelfun
        kernelname
        sample_f
        dsample_f_dx
    end
    methods
        function obj = GPnd(D, theta, kernelname, seed)
            if nargin<4
                seed = 1;
            end
            if nargin <3
                 kernelname = 'Matern52';
            end
            if nargin <2
                theta.cov= [-1;1];
            end
            if nargin <1
                D = 1;
            end
            rng(seed)
            obj.theta = theta;
            obj.D = D;
            obj.xbounds = repmat([0, 1], D, 1);
           
            if strcmp(kernelname, 'Matern52') || strcmp(kernelname, 'Matern32') 
               approximation.method = 'RRGP';
            else
                approximation.method = 'SSGP';
            end
            approximation.nfeatures = 4096;
            approximation.decoupled_bases = 1;
            kernelfun = str2func([kernelname, '_kernelfun']);
            obj.kernelfun = kernelfun;
            obj.kernelname=  kernelname;
            model.kernelname = kernelname;
            model.regularization = 'nugget';
            model.D = D;
            model.kernelfun = kernelfun;
            [obj.sample_f, obj.dsample_f_dx] = sample_GP(obj.theta,  zeros(D,1), [], model, approximation);
        end
        function [y, dydx] = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
           
            y = obj.sample_f(xx);
            
%             if nargout>1
                dydx = obj.dsample_f_dx(xx);
%             end
            
            if strcmp(obj.opt, 'max')
                y = -y;
%                 if nargout>1
                    dydx = -dydx;
%                 end
            end            
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
        end
    end
end


