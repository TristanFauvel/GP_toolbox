function [C, dC, dC_dx] = poly_kernelfun(theta, x0, x, ~, regularization)
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
a =  theta(1);
b =  theta(2);
c =  theta(3);

y0 = x0(2,:);
x0 = x0(1,:);
y = x(2,:);
x = x(1,:);

phi0 = a*x0.*y0+b*x0+c*y0;
phi = a*x.*y+b*x+c*y;

% covariance
C =  (phi0'*phi);
if isequal(x0,x)
    C = (C+C')/2; %to ensure symmetry;
end

%% compute derivative with respect to theta
if nargout>1    %&& nargout<3    
    dC = zeros(n0,n, numel(theta));
    dC(:, :, 1)= 2*a*(x0'.*y0')*(x.*y) + b*x0'*(x.*y) + b*(x0'.*y0')*x + c*y0'*(x.*y)+ c*(x0'.*y0')'*y;
    dC(:, :, 2)= a*x0'*(x.*y)+a*(x0'.*y0')*x+2*b*x0'*x+c*y0'*x+c*x0'*y;
    dC(:, :, 3)= a*(x0'.*y0')*y+a*y0'*(x*y)+b*y0'*x+b*x0'*y+2*c*y0'*y;
end

% if nargout>2
%     
%     dC_dx = zeros(n0,n,n,nd);
%     for i =1:n0
%         for j= 1:n
%             for d= 1:nd
%                 dC_dx(i,j,j,d) = lambda(d).*(x0(d,i)-x(d,j))*C(i,j);
%             end
%         end
%     end
%     
% 
% end

return
