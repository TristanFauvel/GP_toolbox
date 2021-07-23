classdef problem2
    properties
        D = 1;
        xbounds = [2.7,7.5]
        name = 'Problem 2';
        opt = 'max';
        
    end
    methods
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            
            y = sin(x) + sin(10/3*x);
            
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
            if strcmp(obj.opt, 'max')
                y = -y;
            end
        end
    end
end
