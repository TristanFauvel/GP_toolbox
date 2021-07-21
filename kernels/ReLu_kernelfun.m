function [C, dC_dtheta] = ReLu_kernelfun(theta, x0, x, training, regularization)

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
[nd,n]= size(x);
n0 = size(x0,2);
c = theta(1:nd); %c is the inflexion point, where all sample go through
x = x-c;
x0 = x0-c;

sigma2 = theta(end)^2;

arg=x0'*x./(vecnorm(x0,2,1)'*vecnorm(x,2,1));
arg(arg>1)=1;
arg(-arg>1)=-1;
a=acos(arg);

J1=sin(a)+(pi-a).*cos(a);

C= 1/pi*(vecnorm(x0,2,1)'*vecnorm(x,2,1)).*J1;

if training
    C= C+ sigma2*eye(n0);
end
if isequal(x0,x)
    C = (C+C')/2; %to ensure symmetry;
end

% if nargout>1
%     dC_dtheta = zeros(n0,n, numel(theta));
%     dx_c_dtheta = zeros(nd,n,nd);
%     for d =1:nd
%         dx_c_dtheta(d,:,d) = -1;
%     end
%     dC_dx_c = zeros(n0,n,n,nd);
%     dx_c_dxd_c = (x-c)./vecnorm(x-c,2,1);
%     
%     %dadx = -1./sqrt(1+a^2).*
%     for i =1:n0
%         for j= 1:n
%             for d= 1:nd
%                 b = a(i,j);
%                 dx_c_dxd_c = (x(d,j)-c(d))./vecnorm(x(:,j)-c);
%                 dbdxd = (-1/sqrt(1-arg(i,j)^2))*sign(x0(d,i)-c(d))*dx_c_dxd_c;
%                 dC_dx_c(i,j,j,d) = vecnorm(x0(:,i))*(x(d,j)/vecnorm(x(:,j)))*(sin(b)+(pi-b)*sin(b)) + vecnorm(x0(:,i)-c)*vecnorm(x(:,j)-c)*(cos(b)+(pi-b)*cos(b)-sin(b))*dbdxd;
%                 dC_dx_c(i,j,j,d) =  dC_dx_c(i,j,j,d)/pi;
%                 
%                 dC_dtheta(i,j,d) = -2*dC_dx_c(i,j,j,d)*dx_c_dtheta(d,j,d);
%             end
%         end
%     end    
%     dC_dtheta(:,:,end) = 2*theta(end)*eye(n0);
% end
% 
% % dC_dx = zeros(n0,n,n,nd);
% %     %dx_c_dxd = (x(d,j)-c(d))./vecnorm(x(:,i)-c)
% %     dx_c_dxd = (x-c)./vecnorm(x-c,2,1);
% %
% %     %dadx = -1./sqrt(1+a^2).*
% %     for i =1:n0
% %         for j= 1:n
% %             for d= 1:nd
% %                 b = a(i,j);
% %                 dx_c_dxd = (x(d,j)-c(d))./vecnorm(x(:,j)-c);
% %                 dbdxd = (-1/sqrt(1-arg(i,j)^2))*sign(x0(d,i)-c(d))*dx_c_dxd;
% %                 dC_dx(i,j,j,d) = vecnorm(x0(:,i))*(x(d,j)/vecnorm(x(:,j)))*(sin(b)+(pi-b)*sin(b)) + vecnorm(x0(:,i)-c)*vecnorm(x(:,j)-c)*(cos(b)+(pi-b)*cos(b)-sin(b))*dbdxd;
% %             end
% %         end
% %     end
% %     dC_dx =  dC_dx/pi;