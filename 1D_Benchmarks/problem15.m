classdef problem15
    properties
        D = 1;
xbounds = [-5,5];
        name = 'Problem 15';
        opt = 'max';
        
    end
    methods
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            y=(xx.^2-5*xx+6)./(xx.^2+1);            
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
            if strcmp(obj.opt, 'max')
                y = -y;
            end
        end
    end
end
