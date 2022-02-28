classdef Tesch2
    properties
        D = 1
        xbounds = [0 10];
        name = 'Tesch2';
        opt = 'max';
        mean
        var
        takelog
        rescaling  

    end
    methods
        function obj = Tesch2()
            obj.rescaling = 0;             
        end
        
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            x = xx(:);
           
            y = sin(5*x/4) - cos(3*(x-1)/5)/20 - (5*x.^3 + 54*x.^2 - 179*x +159)./100;
            y = norminv(y);

            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
            
           if obj.rescaling
                if obj.takelog
                    if any(y<0)
                        error('Log of negative value')
                    end
                    y = log(y);
                end
                y = (y- obj.mean)./sqrt(obj.var);
            end
            
            if strcmp(obj.opt, 'min')
                y = -y;
            end
        end
    end
end
