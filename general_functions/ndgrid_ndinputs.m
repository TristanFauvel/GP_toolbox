function output = ndgrid_ndinputs(x)
%x = [3,3,4];
%d = 3;

ndim = numel(x);
I=cell(ndim,1);
% construct the neighborhood
for di=1:ndim
    I{di}=1:x(di);
end
[I{1:ndim}]=ndgrid(I{:});
output = [];
for di = 1:ndim
    m = I{di};
    output = [output,m(:)];
end
