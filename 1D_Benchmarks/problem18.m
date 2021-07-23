classdef problem18
    properties
        D = 1;
        xbounds= [0,6];
        name = 'Problem 18';
        opt = 'max';
        
    end
    methods
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            if x<=3
                y=(x-2).^2;
            else
                y=2*log10(x-2)+1;
            end
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
            if strcmp(obj.opt, 'max')
                y = -y;
            end
        end
    end
end
