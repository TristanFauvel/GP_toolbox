function [best, minval] = multistart_minConf(fun, lb, ub, ncandidates, init_guess, options)

%% Multistart seach with minFunc
DEFAULT('ncandidates',10);
DEFAULT('options', []);
best= [];
nd = numel(lb);
starting_points = lb(:)+(ub(:)-lb(:)).*rand(nd,ncandidates);

if ~isempty(init_guess)
    starting_points(:,1)= init_guess;
end
%%
minval=inf;

options.verbose = 0;
x = [];
for k = 1:ncandidates    
    try
        x = minConf_TMP(@(x)fun(x), starting_points(:,k), lb(:), ub(:), options);
        if any(isnan(x(:)))
            error('x is NaN')
        end
        
        val = fun(x);
        if val < minval
            best = x;
            minval = val;
        end
    catch
        disp('Failure')
    end
end

if isempty(x)
    error('Optimization failed')
end
return

