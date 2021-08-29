%From http://infinity77.net/global_optimization/test_functions_1d.html
classdef problem14
    properties
        D = 1;
        xbounds = [0,4];
        name = 'Problem 14';
        opt = 'max';
        
    end
    methods
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            y=-exp(-xx).*sin(2*pi*xx);
            
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
            if strcmp(obj.opt, 'max')
                y = -y;
            end
        end
    end
end
