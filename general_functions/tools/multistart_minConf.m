function [best_x, best_val] = multistart_minConf(fun, lb, ub, ncandidates, init_guess, options, varargin)
opts = namevaluepairtostruct(struct( ...
    'dims', [], ...
    'objective', 'min' ...
    ), varargin);

UNPACK_STRUCT(opts, false)

%% Multistart seach with minFunc
DEFAULT('ncandidates',10);
DEFAULT('options', []);
best_x= [];

x0 = lb;
f = @(x) subfun(x, fun,dims, x0, objective);
if ~isempty(dims)
    lb = lb(dims);
    ub = ub(dims);
end

starting_points = rand_interval(lb, ub, 'nsamples', ncandidates);
if ~isempty(init_guess)
    starting_points(:,1)= init_guess;
end

if strcmp(objective, 'max')
    best_val=-inf;
elseif strcmp(objective, 'min')
    best_val=inf;
end

options.verbose = 0;
x = [];
for k = 1:ncandidates
    try
        x = minConf_TMP(@(x)f(x), starting_points(:,k), lb(:), ub(:), options);
        if any(isnan(x(:)))
            error('x is NaN')
        end
        val = f(x);

        if strcmp(objective, 'max')
             val = -val;
        end

        if (strcmp(objective, 'min') && val < best_val) || (strcmp(objective, 'max') && best_val < val)
            best_x = x;
            best_val = val;
        end
    catch e %e is an MException struct        
        fprintf(1,e.message);
    end
end

if isempty(x)
    error('Optimization failed')
end

if ~isempty(dims)
x0(dims) = best_x;
best_x = x0;
end
end

function [f, df] = subfun(x, fun,dims, x0, objective)

if ~isempty(dims)
    x0(dims) = x;
    [f, df] = fun(x0);
    df = df(dims,:);
else
    [f, df] = fun(x);
end
df= df(:);

if strcmp(objective, 'max')
    f = -f;
    df = -df;
end
end
