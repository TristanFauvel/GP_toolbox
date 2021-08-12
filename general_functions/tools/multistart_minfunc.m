function best = multistart_minfunc(fun, lb, ub, ncandidates, options)

%% Multistart seach with minFunc
options=[];
DEFAULT('ncandidates',10);
DEFAULT('options', []);
best= [];
nD = model.D;
starting_points= NaN(nd,ncandidates);
for i = 1:nd
    starting_points(i,:)= lb(i)+(ub(i)-lb(i))*rand(1,ncandidates);
end

minval=inf;
for k = 1:ncandidates    
    x = minFunc(@(x)fun(x), starting_points(:,k), options);
    val = fun(x);
    if val < minval
        best = x;
        minval = val;
    end
end
return

