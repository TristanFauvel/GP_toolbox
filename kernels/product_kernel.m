function [C, dC, dC_dx, dCxx_dx] = product_kernel(theta, x0, x, training, regularization, kernelfun1, kernelfun2, d1, d2, dtheta1, dtheta2)
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
% dkxx_dx dervative of k(x,x) w.r.t x
%DEFAULT('x', x0);
x01 = x0(d1,:);
x02 = x0(d2,:);
x1 = x(d1,:);
x2 = x(d2,:);
theta1 = theta(dtheta1);
theta2 = theta(dtheta2);

ntr = size(x0,2);
ntst = size(x,2);
% [C, dC, dC_dx] = sum_kernel(theta1, theta2, x01, x1, x02, x2, kernelfun1, kernelfun2, training, regularization)
% DEFAULT('x1', x01);
% DEFAULT('x2', x02);


if nargout == 1
    C1 = kernelfun1(theta1, x01, x1, training, regularization);
    C2 = kernelfun2(theta2, x02, x2, training, regularization);
elseif nargout == 2
    [C1, dC1] = kernelfun1(theta1, x01, x1, training, regularization);
    [C2, dC2] = kernelfun2(theta2, x02, x2, training, regularization);  
    dC =  NaN(ntr, ntst, numel(theta));
    dC(:,:,dtheta1) = dC1.*C2;
    dC(:,:,dtheta2) = dC2.*C1; %note: here we assume that the hyperparameters are independent between the two kernels, which may not always be the case
else 
    [C1, dC1, dC1_dx] = kernelfun1(theta1, x01, x1, training, regularization);
    [C2, dC2, dC2_dx] = kernelfun2(theta2, x02, x2, training, regularization);    
    dC =  NaN(ntr, ntst, numel(theta));
    dC(:,:,dtheta1) = dC1.*C2;
    dC(:,:,dtheta2) = dC2.*C1; %note: here we assume that the hyperparameters are independent between the two kernels, which may not always be the case
    dC_dx = NaN(ntr,ntst,ntst, numel(union(d1,d2)));
    dC_dx(:,:,:,d1) = dC1_dx.*C2;
    dC_dx(:,:,:,d2) = dC2_dx.*C1; %note: here we assume that the kernels act on different dimensions
    %dC_dx(:,:,:,intersect(d1,d2)) = dC1_dx
end

if nargout>3
    dCxx_dx = zeros(1,size(x0,1)); %Only true for stationary covariance functions
end
C = C1.*C2;
if isequal(x0,x)
    C = (C+C')/2; %to ensure symmetry;
end

return