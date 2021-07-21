function [C, dC, dC_dx] = polynomial_kernelfun(theta, x0, x, training, varargin)
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

if numel(theta) ~= 2
    error('The polynomial kernel requires 2 hyperparameters')
end

DEFAULT('regularization', 'nugget'); % do not regularize the base kernels


if nargin==2
    x = x0;
end

[nd, n0] = size(x0);
n = size(x, 2);

deg = 2; %theta(3);
k0 = exp(theta(1));
c = exp(theta(2));
C0 = (x0'*x+c).^deg;
C  = k0*C0;
if isequal(x0,x)
    C = (C+C')/2; %to ensure symmetry;
end

if strcmp(regularization, 'nugget')
    C= nugget_regularization(C);
end
%% compute derivative
if nargout>1   % && nargout<3
    dC = zeros(n0,n, 2);

    dC(:, :, 1) = k0*C0;

    dC(:, :, 2) =   c*k0*deg*(x0'*x+c).^(deg-1);

%     dC(:, :, 3) =  0;
end



if nargout>2

%     dC_dx = zeros(n0,n,n,nd);
%     for i =1:n0
%         for j= 1:n
%             for d= 1:nd
%             dC_dx(i,j,j,d) = x(d,j)*k0*deg*(x0(d,i)'*x(d,j)+c).^(deg-1);
%             end
%         end
%     end

    arg  = k0*deg*(x0'*x+c).^(deg-1);
    dC_dx = zeros(n0,n,n,nd);
    for i =1:n0
        for j= 1:n
            for d= 1:nd
                dC_dx(i,j,j,d) = x0(d,i)*arg(i,j);

            end
        end
    end
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
