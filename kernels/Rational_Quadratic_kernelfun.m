function [C, dC, dC_dx] = Rational_Quadratic_kernelfun(theta, x0, x, training, regularization)
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

if nargin==2
    x = x0;
end

% unpack hyperparameters
[nd, n0] = size(x0);
n = size(x, 2);
if numel(theta)~=3
    error("The rational quadratic kernels requires 3 hyperparameters")
end
lambda       =  exp(theta(1));
k0           =  exp(theta(2));
alpha = exp(theta(3));

dx2 =pdist2(x0',x').^2;
dx = 0.5*lambda*dx2;

C0 =(1+dx/alpha).^(-alpha);
if strcmp(regularization, 'nugget')
    C0= nugget_regularization(C0);
end

% covariance
C =   k0*C0;
if isequal(x0,x)
    C = (C+C')/2; %to ensure symmetry;
end

%% compute derivative
if nargout>1    %&& nargout<3
    dC = zeros(n0,n, numel(theta));
    dC(:, :, 1) = -k0*((1+dx/alpha).^(-alpha-1)).*dx;
    dC(:, :, 2) = C;
    dC(:, :, 3) = alpha*C.*(-log(1+dx/alpha)+dx./(alpha+dx));

end

if nargout>2

    dC_dx = zeros(n0,n,nd);
    if ~isequal(x0,x)
        M =  lambda*C./(1+dx/alpha);
        for i =1:n0
            for j= 1:n
                for d= 1:nd
                    dC_dx(i,j,d) = M(i,j).*(x0(d,i)-x(d,j));
                end
            end
        end
     end
end
return
