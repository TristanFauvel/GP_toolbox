%From http://infinity77.net/global_optimization/test_functions_1d.html
classdef problem13
    properties
        D = 1;
        xbounds = [0.001, 0.99];
        name = 'Problem 13';
        opt = 'max';
        
    end
    methods
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            y=-xx.^(2/3)-(1-xx.^2).^(1/3);
            
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
            if strcmp(obj.opt, 'max')
                y = -y;
            end
        end
    end
end
