function [C, dC, dC_dx] = linear_kernelfun_wnoise(theta, x0, x, training, regularization)
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

if numel(theta) ~= 4
    error('The linear kernel with noise requires 4 hyperparameters')
end

if nargin==2
    x = x0;
end
    
[nd, n0] = size(x0);
n = size(x, 2);



C0 = (x0- theta(2))'*(x- theta(2));
C  = exp(theta(1))*C0+exp(theta(3));

if training %%Corresponds to the situation in which we need to compute the kernel for the training set, for which we take into account the measurement noise.
    C = C + exp(theta(4))*eye(n0);
end

if isequal(x0,x)
    C = (C+C')/2; %to ensure symmetry;
end


%% compute derivative
if nargout>1    && nargout<3
    dC = zeros(n0,n, 2);
    
    dC(:, :, 1) = exp(theta(1))*C0;

    dC(:, :, 2) =   exp(theta(1))*(theta(2)-x-x0);
    
    dC(:, :, 3) =   exp(theta(3));
    
    dC(:, :, 4) =   exp(theta(4))*eye(n0);
end



if nargout>2
    dC = 0;

    xtemp0 = permute((x0- theta(2)), [2 3 1]);
    xtemp  =  permute((x- theta(2)), [2 3 1]);
    dx = (xtemp0-permute(xtemp, [2 1, 3]));
    dC_dx = exp(theta(1))*dx;
    
end


% C0 = x0'*x;
% C = theta(1)^2*C0 + theta(2)^2;
% 
% 
% %% compute derivative
% if nargout>1    && nargout<3
%     dC = zeros(n0,n, 2);
%     
%     dC(:, :, 1) = 2*theta(1)*C0;
% 
%     dC(:, :, 2) =   2*theta(2);
%     
%    
% end
% 
% 
% 
% if nargout>2
%     dC = 0;
% 
%     xtemp0 = permute(x0, [2 3 1]);
%     xtemp  =  permute(x, [2 3 1]);
%     dx = (xtemp0-permute(xtemp, [2 1, 3]));
%     dC_dx = theta(1)^2*dx;
%     
% end

return
