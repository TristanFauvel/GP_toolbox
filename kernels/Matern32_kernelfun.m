function [C, dC, dC_dx] = Matern32_kernelfun(theta, x0, x, ~, regularization)
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

if numel(theta) ~=2
    error('The Matérn 3/2 kernel requires 2 hyperparameters')
end

rho  =  exp(theta(1));
k0           =  exp(theta(2));

% compute sum_i nu_i^2(x_i - x_i')^2
x02 = sum(x0.*x0, 1);
x2  = sum(x.*x, 1);

% dx2 = (x02'+ x2 - 2*x0'*x);
% dx2(dx2<0) = 0;
% r = sqrt(dx2);

r = pdist2(x0',x');
% nu = 5/2;

% C0 = (1 + sqrt(5)*r/l + 5*r.^2/(3*l^2)).*exp(-sqrt(5)*r/l);
C0 = (1 + sqrt(3)*r/rho).*exp(-sqrt(3)*r/rho);
if strcmp(regularization, 'nugget')
    C0= nugget_regularization(C0);
end

% covariance
C =   k0*C0;
if isequal(x0,x)
    C = (C+C')/2; %to ensure symmetry;
end

%% compute derivative
if nargout>1   % && nargout<3
    dC = zeros(n0,n, numel(theta));
    dC(:, :, 1) = k0*exp(-sqrt(3)*r/rho).*(3*r.^2/rho^2);  % k0*(sqrt(3)*r/rho^2*C + (1-sqrt(5)*r/rho^2 - 10/2*r^2/rho^3)*exp(-sqrt(3)*r/rho));%-(nu+1)*C/rho;
    dC(:, :, 2) = k0*C0;
end

if nargout>2
    
    %      xtemp0 = permute(x0, [2 3 1]);
    %     xtemp  =  permute(x, [2 3 1]);
    %     dx = (xtemp0-permute(xtemp, [2 1, 3]));
    %
    dC_dr = sqrt(3)/rho*k0*exp(-sqrt(3)*r/rho)-sqrt(3)/rho*C;
    %      dC_dx = dC_dr*dx_dr;
    
    %     xtemp0 = permute(x0, [2 3 1]);
    %     xtemp  =  permute(x, [2 3 1]);
    %     dx = (xtemp0-permute(xtemp, [2 1, 3]));
    %     dr_dx = - sign(x0 - x);
    %     dC_dx =0;% dC_dr*dr_dx;
    %dC = 0;
    
    %
    %     xtemp0 = permute(x0, [2 3 1]);
    %     xtemp  =  permute(x, [2 3 1]);
    %     dx = (xtemp0-permute(xtemp, [2 1, 3]));
    %     dC_dx = k0*(sqrt(5)/rho +10*r/(3*rho^2)).*exp(-sqrt(5)*r/rho) - sqrt(5)/rho*C; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FAUX
    dC_dx = zeros(n0,n,n,nd);
    
    if ~isequal(x0,x)
        if nd>1
            for i =1:n0
                for j= 1:n
                    if r(i,j)~=0
                        for d= 1:nd
                            dC_dx(i,j,j,d) = dC_dr(i,j)*(x(d,j)-x0(d,i))./(r(i,j)); % lambda(d).*(x0(d,i)-x(d,j))*C(i,j);
                        end
                    end
                end
            end
        else
            for i =1:n0
                for j= 1:n
                    if r(i,j)~=0
                        dC_dx(i,j,j) = dC_dr(i,j)*(x(j)-x0(i))./(r(i,j)); % lambda(d).*(x0(d,i)-x(d,j))*C(i,j);
                    end
                end
            end
        end
    end
end
return
