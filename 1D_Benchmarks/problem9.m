classdef problem9
    properties
        D = 1;
        xbounds = [3.1,20.4];
        name = 'Problem 9';
        opt = 'max';
        
    end
    methods
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            y =sin(xx) + sin(2/3*xx);
            
            
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
            if strcmp(obj.opt, 'max')
                y = -y;
            end
        end
    end
end