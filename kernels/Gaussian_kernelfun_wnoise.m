function [C, dC, dC_dx] = Gaussian_kernelfun_wnoise(theta, x0, x,training, regularization)
%[dx2, dC, dC_dx] = Gaussian_kernelfun_wnoise(theta, x0, x)

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
% x    [nd, N_tst]:  training data
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
lambda       =  exp(theta(1));
k0           =  exp(theta(2));
sigma2       =  exp(theta(3));

% compute sum_i nu_i^2(x_i - x_i')^2
x02 = sum(x0.*x0, 1);
x2  = sum(x.*x, 1);

dx2 = (x02'/2 + x2/2 - x0'*x);

% take exponential
C0 = exp( -lambda*dx2 );


% % covariance
% if nargin==2 || isequal(x0,x) %Corresponds to the situation in which we need to compute the kernel for the training set.
%     C =   k0*C0 + sigma2*eye(n0);
% elseif nargin == 3
%     C =   k0*C0;
% end

if training %%Corresponds to the situation in which we need to compute the kernel for the training set, for which we take into account the measurement noise.
    C =   k0*C0 + sigma2*eye(n0);
else %if ~training && isequal(x0,x)
    C =   k0*C0;
end

if strcmp(regularization, 'nugget')
    C= nugget_regularization(C);
end

if isequal(x0,x)
    C = (C+C')/2; %to ensure symmetry;
end

%% compute derivative
if nargout>1    && nargout<3
    
    dC = zeros(n0,n, 2);
    
    dC(:, :, 1) = -lambda.*dx2.*k0.*C0;
    
    dC(:, :, 2) =   k0*C0;
    
    %     if isequal(x0,x)
    %         dC(:, :, 3) =   sigma2*eye(n0);
    %     end
    if training
        dC(:, :, 3) =   sigma2*eye(n0);
    end
    
end



if nargout>2
     dC_dx = zeros(n0,n,n,nd);
         if ~isequal(x0,x)       
    for i =1:n0
        for j= 1:n
            dC_dx(i,j,j,:) = lambda.*(x0(:,i)-x(:,j))*C(i,j);
        end
    end
         end
%     
%     xtemp0 = permute(x0, [2 3 1]);
%     xtemp  =  permute(x, [2 3 1]);
%     dx = (xtemp0-permute(xtemp, [2 1, 3]));
%     
%     dC_dx = lambda*k0*dx.*C0;
    
    dC = 0;
    
end

return
