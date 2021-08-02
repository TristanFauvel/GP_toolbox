function [C, dC, dC_dx] = sum_kernel(theta, x0, x, training, kernelfun1, kernelfun2, d1, d2, dtheta1, dtheta2, regularization)
DEFAULT('x', x0);
x01 = x0(d1,:);
x02 = x0(d2,:);
x1 = x(d1,:);
x2 = x(d2,:);
theta1 = theta(dtheta1);
theta2 = theta(dtheta2);

ntr = size(x0,2);
ntst = size(x,2);
% [C, dC, dC_dx] = sum_kernel(theta1, theta2, x01, x1, x02, x2, kernelfun1, kernelfun2, training)
% DEFAULT('x1', x01);
% DEFAULT('x2', x02);

dC1_dx = [];
dC2_dx = [];

reg= 'none';
C1 = kernelfun1(theta1, x01, x1, training, reg);
    C2 = kernelfun2(theta2, x02, x2, training, reg);
if nargout == 2
    [C1, dC1] = kernelfun1(theta1, x01, x1, training, reg);
    [C2, dC2] = kernelfun2(theta2, x02, x2, training, reg);  
    dC =  NaN(ntr, ntst, numel(theta));
    dC(:,:,dtheta1) = dC1;
    dC(:,:,dtheta2) = dC2; %note: here we assume that the hyperparameters are independent between the two kernels, which may not always be the case
elseif nargout == 3
    [C1, dC1, dC1_dx] = kernelfun1(theta1, x01, x1, training, reg);
    [C2, dC2, dC2_dx] = kernelfun2(theta2, x02, x2, training, reg);    
    dC =  NaN(ntr, ntst, numel(theta));
    dC(:,:,dtheta1) = dC1;
    dC(:,:,dtheta2) = dC2; %note: here we assume that the hyperparameters are independent between the two kernels, which may not always be the case
    dC_dx = NaN(ntr,ntst,ntst, numel(union(d1,d2)));
    dC_dx(:,:,:,d1) = dC1_dx;
    dC_dx(:,:,:,d2) = dC2_dx; %note: here we assume that the kernels act on different dimensions
    %dC_dx(:,:,:,intersect(d1,d2)) = dC1_dx
end
C = C1 + C2;

if strcmp(regularization, 'nugget')
    C = nugget_regularization(C);
end

return