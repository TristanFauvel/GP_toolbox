classdef Tesch1
    properties
        D = 1
        xbounds = [0 10];
        name = 'Tesch1';
        opt = 'max';
        mean
        var
        takelog
        rescaling  

    end
    methods
        function obj = Tesch1()
            obj.rescaling = 0;             
        end
        
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            x = xx(:);
           
            y = sin(x) - cos(3*x)/4 + (x.^3 - 13*x.^2 - 29*x -55)./50;
            
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
