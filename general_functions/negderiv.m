function [fx, dfdx] = negderiv(x,f,df)
if any(isnan(x))
    warning('x is NaN')
end
% Function that groups f and df to use minFunc
% disp(x)
fx = -f(x);
% if nargout > 1
dfdx = - df(x);
% end