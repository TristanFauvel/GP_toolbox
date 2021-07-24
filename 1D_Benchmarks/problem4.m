classdef problem4
    properties
        D = 1;
        xbounds = [1.9,3.9];
        name = 'Problem 4';
        opt = 'max';
        
    end
    methods
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            y=-(16*xx.^2-24*xx+5).*exp(-xx);
            
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
            if strcmp(obj.opt, 'max')
                y = -y;
            end
        end
    end
end
