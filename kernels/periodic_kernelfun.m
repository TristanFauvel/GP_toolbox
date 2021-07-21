function C = periodic_kernelfun(theta, x0, x, training, regularization)
%% C = covfun(theta, x0)
% compute covariance of outputs
% 
% INPUTS: 
% theta [nd, 1]: hyperparameters
%               theta(1:nd) = nu, for each dimension
%               theta(nd+1) = scaling factor: gaussian kernel
%               theta(nd+2) = ouput noise
% 
% x0   [nd, N_tr]:  training data
% x0   [nd, N_tst]:  training data
%
% OUTPUT
% C  = [N, N]:           covariance of p(y|x)
% dC = [N, N, ntheta]:   derivative of C w.r.t theta 
DEFAULT('regularization', 'nugget'); % do not regularize the base kernels

% opts = namevaluepairtostruct(struct( ...
%     'regularization', 'nugget' ...
%     ), varargin);
% UNPACK_STRUCT(opts, false)
% 
if nargin==2
    x = x0;
end
    
% unpack hyperparameters
[nd, n0] = size(x0);
n = size(x, 2);

if numel(theta) ~= nd+2
    error('The number of hyperparameters for the periodic kernel is nd + 1')
end

lambda       =  exp(theta(1:nd));
k0           =  exp(theta(nd+1));
p = theta(nd+2);
% compute sum_i nu_i^2(x_i - x_i')^2
dx = zeros(n0,n);
for i =1:nd
    dx=dx + 2*lambda(i)*(sin(pi*(abs(x(i,:)-x0(i,:)'))/p)).^2;
end

C0 =exp(-dx);
if strcmp(regularization, 'nugget')
    C0= nugget_regularization(C0);
end

% covariance
C =   k0*C0;
if isequal(x0,x)
    C = (C+C')/2; %to ensure symmetry;
end
return
