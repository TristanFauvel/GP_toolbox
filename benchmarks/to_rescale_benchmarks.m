close all
add_gp_module

seed=1;
maxiter= 1000;

% filename = [pathname, '/1D_Benchmarks/1D_benchmarks_table.mat'];
folder = [pathname, '/Benchmarks/'];
filename = [folder,'benchmarks_table.mat'];

load(filename, 'benchmarks_table')
objectives = benchmarks_table.fName;

N= numel(objectives);

T = benchmarks_table;

nsamps=  50000;

m = zeros(1,N);
v = zeros(1,N);
take_log = zeros(1,N);
for j = 1:N
    rng(seed)
    
    objective = char(objectives(j));
    
    disp(['Function : ' , num2str(j)])
    [g, theta.cov, model] = load_benchmarks(objective, [], benchmarks_table,0);
    
    lb = model.lb;
    ub = model.ub;
    lb_norm = model.lb_norm;
    ub_norm = model.ub_norm;
    hyp_lb = model.hyp_lb;
    hyp_ub = model.hyp_ub;
    
    xsamples = rand_interval(lb,ub,'nsamples', nsamps);
    ysamples= -g(xsamples);
    m(j) = mean(ysamples);
    v(j) = var(ysamples);
    
    if v(j)>1e2 && all(ysamples>0)
        take_log(j) = 1;
        ysamples= log(-g(xsamples));
        m(j) = mean(ysamples);
        v(j) = var(ysamples);
    end
end
Mean = m(:);
Variance = v(:);
TakeLog = take_log(:);
Names = T.Name;
t = table(Names, Mean, Variance, TakeLog);

% save([folder,'benchmarks_rescaling'], 't')
