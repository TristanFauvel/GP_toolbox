clear all
add_gp_module;
fName = {'ackley'; 'beale'; 'boha1'; 'camel3'; 'camel6'; 'colville'; 'crossit'; ...
    'dixonpr'; 'drop'; 'egg'; 'forretal08'; 'goldpr'; 'griewank'; 'grlee12';  ...
    'hart3'; 'hart4'; 'hart6'; 'holder'; 'langer';'levy';'levy13';'perm0db';'permdb';...
    'powell'; 'rosen'; 'rothyp'; 'schaffer4'; 'schwef'; 'shekel'; 'shubert'; ...
    'spheref'; 'sumsqu'; 'trid'; 'Ursem_waves'};

% fName = {'forretal08';'grlee12'; 'problem2'; 'problem3';'problem4';'problem5'; ...
%     'problem6';'problem7';'problem8';'problem9';'problem10';'problem11';'problem12';...
%     'problem13';'problem14';'problem15';'problem18';'problem20';'problem21';};

N = numel(fName);

for i = 1:N
    fun = str2func(fName{i});
    fun = fun(0);
    D(i,1) = fun.D;
    Name{i,1} = fun.name;
    ARD{i,1} = NaN(D(i,1)+1,1);
    Matern32{i,1} = NaN(2,1);
    Matern52{i,1} = NaN(2,1);    
end
benchmarks_table = table(Name, fName, D, ARD, Matern32, Matern52);

% table2latex(benchmarks_table, [pathname, '/benchmarks.tex'])

fName = categorical(fName);
Name = categorical(Name);
benchmarks_table = table(Name, fName, D, ARD, Matern32, Matern52);

Kernel_name = categorical(benchmarks_table.Kernel_name);
benchmarks_table.Kernel_name  = Kernel_name ;

Kernel= categorical(benchmarks_table.Kernel);
benchmarks_table.Kernel  = Kernel ;

save([pathname, '/Benchmarks/benchmarks_table.mat'],'benchmarks_table')
% save([pathname, '/1D_Benchmarks/1D_benchmarks_table.mat'],'benchmarks_table')

% load([pathname, '/Benchmarks/benchmarks_table.mat'],'benchmarks_table')

load([pathname, '/Benchmarks/benchmarks_table.mat'],'benchmarks_table')
% 
% %% To add the bounds to the table
% N = size(benchmarks_table,1);
% tb = zeros(
% for i = 1:N
%      fun = str2func(fName{i});
%      obj = fun(0);
%      bounds= obj.xbounds;
% end
% 
T = benchmarks_table(:,[1,3,8]);
T.Name= char(T.Name);
T.Kernel_name= char(T.Kernel_name);
% table2latex(T, [pathname, '/benchmarks.tex'])





