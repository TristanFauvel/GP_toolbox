function [C, dC, dC_dx] = ReLu_wbias_kernelfun(bias, x0, x)
%[C, dC, dC_dx] = kernelfun(theta, x0, x)
%% C = covfun(theta, x0)
% compute covariance of outputs
%
% INPUTS:
% theta [2, 1]: hyperparameters
%               theta(1) = magnitude of kernel
%               theta(2) = spatial scale
%
% x0   [N, 1]:  training data
% x    [M, 1]:  test data (if empty then x<-x0)
%
% OUTPUT
% C  = [N, M]:           covariance of p(y|x)
% dC = [N, M, ntheta]:   derivative of C w.r.t theta
DEFAULT('regularization', 'nugget'); % do not regularize the base kernels

DEFAULT('x', x0);

C=arc_cosine_kernel(x0,x)+bias;
if isequal(x0,x)
    C = (C+C')/2; %to ensure symmetry;
end

if nargout>1 %&& nargout<3
    dC = 1;
end

if nargout>2
    nd= size(x,1);
    Ntest= size(x,2);
    Ntr=size(x0,2);
    arg=x0'*x./(vecnorm(x0)'*vecnorm(x));
    arg(arg>1)=1;
    arg(-arg>1)=-1;
    theta=acos(arg);
    dC_dx=NaN(Ntr,Ntest, nd);
    for i = 1:Ntest
        for j = 1:Ntr
           dC_dx(j,i,:)= x(:,i)'./vecnorm(x(:,i))^2*C(j,i)+(pi-theta(j,i))./pi.*sign(sin(theta(j,i))).*(x0(:,j)'-x0(:,j)'*x(:,i)*x(:,i)'/(vecnorm(x(:,i))^2));
        end
    end
end
