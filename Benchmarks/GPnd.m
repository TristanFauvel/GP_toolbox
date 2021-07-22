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
                theta= [-1;1];
            end
            if nargin <1
                D = 1;
            end
            rng(seed)
            obj.theta= theta;
            obj.D = D;
            obj.xbounds = repmat([0, 1], D, 1);
           
            if strcmp(kernelname, 'Matern52') || strcmp(kernelname, 'Matern32') 
               approximation_method = 'RRGP';
            else
                approximation_method = 'SSGP';
            end
            nfeatures = 4096;
            decoupled_bases = 1;
            kernelfun = str2func([kernelname, '_kernelfun']);
            obj.kernelfun = kernelfun;
            obj.kernelname=  kernelname;
            obj.sample_f = sample_GP(obj.theta, zeros(D,1), [], kernelname, approximation_method, decoupled_bases,nfeatures, kernelfun);
        end
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            y = obj.sample_f(xx);
            
            if strcmp(obj.opt, 'max')
                y = -y;
            end
            
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
        end
    end
end


