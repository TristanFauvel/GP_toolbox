function best = multistart_minfuncBC(fun, lb, ub, ncandidates, options)

%% Multistart seach with minFunc
DEFAULT('ncandidates',10);
DEFAULT('options', []);
best= [];
nD = model.D;
starting_points= NaN(nd,ncandidates);
% for i = 1:nd
%     range = linspace(lb(i), ub(i), ncandidates+2);
%     starting_points(i,:)= range(2:end-1);
% end

for i = 1:nd
%     range = linspace(lb(i), ub(i), ncandidates+2);
    starting_points(i,:)= lb(i)+(ub(i)-lb(i))*rand(1,ncandidates);
end

%%
minval=inf;
for k = 1:ncandidates    
    x = minFuncBC(@(x)fun(x), starting_points(:,k), lb, ub, options);
    val = fun(x);
    if val < minval
        best = x;
        minval = val;
    end
end
return
