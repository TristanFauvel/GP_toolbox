function [C, dC, dC_dx] = ARD_kernelfun_wnoise(theta, x0, x, training, regularization)
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

if nargin==2
    x = x0;
end

% unpack hyperparameters
[nd, n0] = size(x0);
n = size(x, 2);
lambda       =  exp(theta(1:nd));
k0           =  exp(theta(nd+1));
sigma2       = exp(theta(nd+2));

% compute sum_i nu_i^2(x_i - x_i')^2
dx = zeros(n0,n);
for i =1:nd
    dx=dx + 0.5*lambda(i)*(x(i,:)-x0(i,:)').^2;
end

C0 =exp(-dx);


C =   k0*C0;

if training %%Corresponds to the situation in which we need to compute the kernel for the training set, for which we take into account the measurement noise.
    C =   k0*C0 + sigma2*eye(n0);
else 
    C =   k0*C0;
end
if isequal(x0,x)
    C = (C+C')/2; %to ensure symmetry;
end

%% compute derivative
if nargout>1    && nargout<3
    dC = zeros(n0,n,nd+2);
    
    for i=1:nd
        dx=(x(i,:)-x0(i,:)').^2;
        dC(:, :, i) = -0.5*lambda(i).*dx.*C;
    end
    
    dC(:, :, nd+1) =   k0*C0;
    
    if training
        dC(:, :, nd+2) =   sigma2*eye(n0);
    end
    
end


dC_dx = zeros(n0,n, nd);
dC_dx0 = zeros(n,n0, nd);

if nargout>2
    
       
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
    
%     xtemp0 = permute(x0, [2 3 1]);
%     xtemp  =  permute(x, [2 3 1]);
%     dx = (xtemp0-permute(xtemp, [2 1, 3]));
%     for i =1:nd
%         dC_dx(:,:,i) = lambda(i)*k0*dx(:,:,i).*C0;
%     end
    
%     dx0 = (xtemp-permute(xtemp0, [2 1, 3]));
%     for i =1:nd
%         dC_dx0(:,:,i) = lambda(i)*k0*dx0(:,:,i).*C0;
%     end
    dC = 0; 
end
return
