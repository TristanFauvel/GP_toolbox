function [best, minval] = multistart_minConf(fun, lb, ub, ncandidates, init_guess, options, varargin)
opts = namevaluepairtostruct(struct( ...
    'dims', [] ...
    ), varargin);

UNPACK_STRUCT(opts, false)

%% Multistart seach with minFunc
DEFAULT('ncandidates',10);
DEFAULT('options', []);
best= [];

if ~isempty(dims)
    x0 = lb; 
    fun = @(x) subfun(x, fun,dims, x0);
    lb = lb(dims);
    ub = ub(dims);
end


starting_points = rand_interval(lb, ub, 'nsamples', ncandidates);
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
end

function [f, df] = subfun(x, fun,dims, x0)
    x0(dims) = x;
    [f, df] = fun(x0);    
    df = df(:,dims);
end
