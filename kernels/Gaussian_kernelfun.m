function [C, dC, dC_dx] = Gaussian_kernelfun(theta, x0, x, training, regularization)
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
% x    [nd, N_tst]:  test data
%
% OUTPUT
% C  = [N_tr, N_tst]:           covariance of p(y|x)
% dC = [N_tr, N_tst, ntheta]:   derivative of C w.r.t theta
DEFAULT('regularization', 'nugget')
DEFAULT('x', x0); % do not regularize the base kernels

% opts = namevaluepairtostruct(struct( ...
%     'regularization', 'nugget' ...
%     ), varargin);
% UNPACK_STRUCT(opts, false)
 
 
if isempty(theta)
    C = 2;
    return
end



if any(isnan([x(:);x0(:)]))
    error('x is NaN')
end

if numel(theta) ~= 2
    error('The number of hyperparameters for the Gaussian kernel is 2')
end


% unpack hyperparameters
[nd, n0] = size(x0);
n = size(x, 2);
lambda       =  exp(theta(1));
k0           =  exp(theta(2));

% compute sum_i nu_i^2(x_i - x_i')^2
x02 = sum(x0.*x0, 1);
x2  = sum(x.*x, 1);

dx2 = (x02'/2 + x2/2 - x0'*x);
dx2(0>dx2)=0;

% take exponential
C0 = exp( -lambda*dx2 );

if strcmp(regularization, 'nugget')
    C0= nugget_regularization(C0);
end
% covariance
C =   k0*C0 ;
if isequal(x0,x)
    C = (C+C')/2; %to ensure symmetry;
end

%% compute derivative
if nargout>1 %&& nargout<3
    
    dC = zeros(n0,n, 2);
    
    dC(:, :, 1) = -lambda.*dx2.*C;
    
    dC(:, :, 2) =   k0*C0;
    
end

if nargout>2
    %      xtemp0 = permute(x0, [2 3 1]);
    %     xtemp  =  permute(x, [2 3 1]);
    %     dx = (xtemp0-permute(xtemp, [2 1, 3]));
    %
    %     dC_dx = lambda.*dx.*C;
    %     dC = 0;
%     if training
%         dC_dx = zeros(n,n,n,nd);
%     else

        dC_dx = zeros(n0,n,n,nd);
            if ~isequal(x0,x)       

        for i =1:n0
            for j= 1:n
                dC_dx(i,j,j,:) = lambda.*(x0(:,i)-x(:,j))*C(i,j);
            end
        end
            end
%     end
    
    % dC = 0;
end
return