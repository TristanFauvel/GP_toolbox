function [y, dydx] = logistic(x)
%LOGISTIC computes 1./(1+exp(-x)) elementwise

y = 1./(1+exp(-x));

if nargout>1
    dydx= y.*(1-y);
end