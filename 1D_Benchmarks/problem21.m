classdef problem21
    properties
        D = 1;
        xbounds= [0,10];
        name = 'Problem 21';
        opt = 'max';
        
    end
    methods
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            y = xx.*sin(xx) + xx.*cos(2*xx);
            
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
            if strcmp(obj.opt, 'max')
                y = -y;
            end
        end
    end
end