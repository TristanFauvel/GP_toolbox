function [C, dC, dC_dx] = ARD_kernelfun(theta, x0, x, ~, regularization)
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

if nargin==2
    x = x0;
end
    
% unpack hyperparameters
[nd, n0] = size(x0);
n = size(x, 2);


if isempty(theta)
    C = nd+1;
    return
end
if numel(theta) ~= nd+1
    error('The number of hyperparameters for the ARD kernel is nd + 1')
end

lambda       =  exp(theta(1:nd));
k0           =  exp(theta(nd+1));

% compute sum_i nu_i^2(x_i - x_i')^2
dx = zeros(n0,n);
for i =1:nd
    dx=dx + 0.5*lambda(i)*(x(i,:)-x0(i,:)').^2;
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

%% compute derivative
if nargout>1    
    dC = zeros(n0,n, numel(theta));
    for i=1:nd
        dx=(x(i,:)-x0(i,:)').^2;
        dC(:, :, i) = -0.5*lambda(i).*dx.*C;
    end
    dC(:, :, nd+1) = k0*C0;   
end

if nargout>2
%     dC_dx = zeros(n0,n,nd);
%     
%     if ~isequal(x0,x)
%         for i =1:n0
%             for j= 1:n
%                 for d= 1:nd
%                     dC_dx(i,j,d) = lambda(d).*(x0(d,i)-x(d,j))*C(i,j);
%                 end
%             end
%         end
%     end

    dC_dx = zeros(n0,n,n,nd);
    
    if ~isequal(x0,x)       
    for i =1:n0
        for j= 1:n
            for d= 1:nd
                dC_dx(i,j,j,d) = lambda(d).*(x0(d,i)-x(d,j))*C(i,j);
            end
        end
    end
    end
end

return
