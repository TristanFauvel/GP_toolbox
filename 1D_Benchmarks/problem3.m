%From http://infinity77.net/global_optimization/test_functions_1d.html
classdef problem3
    properties
        D = 1;
        xbounds = [-10,10];
        name = 'Problem 3';
        opt = 'max';
        
    end
    methods
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            
            k = 1:6;
            y = -sum(k'.*sin((k+1)'.*xx+k'),1);
            
            
            
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
            if strcmp(obj.opt, 'max')
                y = -y;
            end
        end
    end
end
