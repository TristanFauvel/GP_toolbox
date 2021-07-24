classdef problem6
    properties
        D = 1;
        xbounds = [-10,10];
        name = 'Problem 6';
        opt = 'max';
        
    end
    methods
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            
            y = -(xx+sin(xx)).*exp(-xx.^2);
            
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
            if strcmp(obj.opt, 'max')
                y = -y;
            end
        end
    end
end
