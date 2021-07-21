function x = logistic(x)
%LOGISTIC computes 1./(1+exp(-x)) elementwise
%
% Maps an unconstrained input to [0,1]

% Iain Murray, November 2007

x = 1./(1+exp(-x));
