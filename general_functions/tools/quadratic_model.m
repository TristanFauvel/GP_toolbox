function [y,dydx] = quadratic_model(x, weight) 

[phi, dphidx]= quadratic_features(x);

y = -weight'*phi;
dydx= mtimesx(weight',dphidx);
dydx= -squeeze(dydx);