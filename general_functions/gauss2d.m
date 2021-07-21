function [val, exponent]= gauss2d(a, x, y, sigma, center)
pos= [x(:), y(:)]';

% exponent = diag((pos-center)'*inv(sigma)*(pos-center)/2);
% exponent=reshape(exponent, size(x));
% val       = a*(exp(-exponent));
val = a*mvnpdf(pos',center',eye(2)*sigma);
val = reshape(val, size(x));
return 


