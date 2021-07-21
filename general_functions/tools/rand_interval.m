function sample = rand_interval(lb,ub, varargin)

opts = namevaluepairtostruct(struct( ...
  'nsamples', 1 ...
    ), varargin);

UNPACK_STRUCT(opts, false)

n = numel(lb);
sample = NaN(n,nsamples);
for i = 1:n
    sample(i,:) = lb(i) + (ub(i)-lb(i))*rand(1,nsamples);
end