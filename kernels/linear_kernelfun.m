function [C, dC, dC_dx] = linear_kernelfun(theta, x0, x, ~, regularization)
%[dx2, dC, dC_dx] = linear_kernelfun(theta, x0, x)

%% C = covfun(theta, x0)
% compute covariance of outputs
% 
% INPUTS: 
% theta= []: hyperparameters
%               
% 
% x0   [nd, N_tr]:  training data
% x0   [nd, N_tst]:  training data
%
% OUTPUT
% C  = [N, N]:           covariance of p(y|x)
% dC = [N, N, ntheta]:   derivative of C w.r.t theta 
DEFAULT('regularization', 'nugget'); % do not regularize the base kernels

if numel(theta) ~= 3
    error('The linear kernel requires 3 hyperparameters')
end

if nargin==2
    x = x0;
end
    
[nd, n0] = size(x0);
n = size(x, 2);


C0 = (x0- theta(2))'*(x- theta(2));
C  = exp(theta(1))*C0+exp(theta(3));
if isequal(x0,x)
    C = (C+C')/2; %to ensure symmetry;
end


%% compute derivative
if nargout>1   % && nargout<3
    dC = zeros(n0,n, 2);
    
    dC(:, :, 1) = exp(theta(1))*C0;

    dC(:, :, 2) =  exp(theta(1))*(-(x' - theta(2))*ones(nd,1) -ones(nd,1)'*(x-theta(2))); % exp(theta(1))*(2*theta(2)-x-x0');
    
    dC(:, :, 3) =   exp(theta(3));
end



if nargout>2      
    dC_dx = zeros(n0,n,n,nd);
    for i =1:n0
        for j= 1:n
              dC_dx(i,j,j,:) =  (x0(:,i)-theta(2));
        end
    end
    dC_dx = exp(theta(1))*dC_dx;
end

if strcmp(regularization ,'nugget')
    C = nugget_regularization(C);
end


return
