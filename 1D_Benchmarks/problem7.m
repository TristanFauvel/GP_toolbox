classdef problem7
    properties
        D = 1;
        xbounds = [2.7,7.5];
        name = 'Problem 7';
        opt = 'max';
        
    end
    methods
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            
            y = sin(xx)+sin(10*xx/3) + log10(xx) -0.84*xx+3;
            
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
            if strcmp(obj.opt, 'max')
                y = -y;
            end
        end
    end
end
